use nalgebra::{DMatrix, DVector, Matrix6, Vector3, Vector6, Unit};
use anyhow::{Result, anyhow, Context};
use rand_distr::{Normal, Poisson, Distribution};
use rand::prelude::*;
use std::fs::File;
use std::io::{Write, BufRead, BufReader};

const C_LIGHT: f64 = 299_792.458;

fn parse_hms(s: &str) -> Result<f64> {
    let parts: Vec<&str> = s.split(':').collect();
    if parts.len() < 3 { return Ok(0.0); }
    let h: f64 = parts[0].parse()?;
    let m: f64 = parts[1].parse()?;
    let s: f64 = parts[2].parse()?;
    Ok((h + m/60.0 + s/3600.0) * 15.0)
}

fn parse_dms(s: &str) -> Result<f64> {
    let parts: Vec<&str> = s.split(':').collect();
    if parts.len() < 3 { return Ok(0.0); }
    let d: f64 = parts[0].parse()?;
    let m: f64 = parts[1].parse()?;
    let s: f64 = parts[2].parse()?;
    let sign = if d < 0.0 || parts[0].starts_with('-') { -1.0 } else { 1.0 };
    Ok(sign * (d.abs() + m/60.0 + s/3600.0))
}

const GM_SUN: f64 = 132_712_440_018.0;

#[derive(Debug, Clone)]
#[allow(dead_code)]
struct Pulsar {
    id: String,
    dir: Unit<Vector3<f64>>,
    flux: f64,
    width: f64,
    period: f64,
}

impl Pulsar {
    fn new(id: &str, ra_deg: f64, dec_deg: f64, period: f64, flux: f64) -> Self {
        let ra = ra_deg.to_radians();
        let dec = dec_deg.to_radians();
        let x = dec.cos() * ra.cos();
        let y = dec.cos() * ra.sin();
        let z = dec.sin();
        Self {
            id: id.to_string(),
            dir: Unit::new_normalize(Vector3::new(x, y, z)),
            flux,
            width: 0.05,
            period,
        }
    }

    fn from_par_file(path: &str) -> Result<Self> {
        let file = File::open(path).with_context(|| format!("Erro ao abrir {}", path))?;
        let reader = BufReader::new(file);

        let mut id = String::new();
        let mut ra = 0.0;
        let mut dec = 0.0;
        let mut period = 0.010;
        
        let mut flux = 1.0; 

        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') || line.trim().is_empty() { continue; }
            
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() < 2 { continue; }

            match parts[0] {
                "PSR" | "PSRJ" => id = parts[1].to_string(),
                "RAJ" => ra = parse_hms(parts[1])?,
                "DECJ" => dec = parse_dms(parts[1])?,
                "F0" => {
                    let f0: f64 = parts[1].parse().unwrap_or(1.0);
                    period = 1.0 / f0;
                },
                "P0" => period = parts[1].parse().unwrap_or(0.010),
                _ => {}
            }
        }

        if id.contains("1937") { flux = 5.0; }
        else if id.contains("0437") { flux = 8.0; }

        Ok(Self::new(&id, ra, dec, period, flux))
    }

    fn simulate_photons(&self, duration: f64, rng: &mut impl rand::Rng) -> Vec<f64> {
        let expected_photons = self.flux * duration;
        if expected_photons <= 0.0 { return Vec::new(); }

        let poisson = Poisson::new(expected_photons).unwrap();
        let n_photons = poisson.sample(rng) as usize;
        
        let mut timestamps = Vec::with_capacity(n_photons);

        let pulse_phase = Normal::new(0.5, self.width / 2.355).unwrap();

        for _ in 0..n_photons {
            let is_signal = rng.gen_bool(0.3); 
            
            let phase = if is_signal {
                let mut ph = pulse_phase.sample(rng);
                ph = ph - ph.floor();
                ph
            } else {
                rng.gen::<f64>()
            };

            let t_window = rng.gen::<f64>() * duration;
            let cycles = (t_window / self.period).floor();
            let final_time = (cycles + phase) * self.period;

            if final_time < duration {
                timestamps.push(final_time);
            }
        }
        
        timestamps.sort_by(|a, b| a.partial_cmp(b).unwrap());
        timestamps
    }

    fn shapiro_delay(&self, pos_spacecraft: &Vector3<f64>) -> f64 {
        let r_vec = pos_spacecraft; 
        let r = r_vec.norm();
        if r < 1.0 { return 0.0; }

        let cos_theta = r_vec.normalize().dot(&self.dir);

        let factor = (2.0 * GM_SUN) / C_LIGHT.powi(3);
        
        -factor * (1.0 - cos_theta).ln()
    }

    #[allow(dead_code)]
    fn calculate_observation_noise(&self, integration_time: f64) -> f64 {
        if integration_time <= 0.0 { return 1.0; }
        let n_ph = self.flux * integration_time;
        
        (self.width * self.period) / (2.0 * n_ph.sqrt())
    }
}

struct KalmanFilter {
    state: Vector6<f64>,
    covariance: Matrix6<f64>,
}

impl KalmanFilter {
    fn new(initial_pos: Vector3<f64>) -> Self {
        let mut state = Vector6::zeros();
        state.fixed_rows_mut::<3>(0).copy_from(&initial_pos);
        
        let covariance = Matrix6::identity() * 1000.0; 

        Self { state, covariance }
    }

    fn predict(&mut self, dt: f64) {
        let mut f = Matrix6::identity();
        f[(0, 3)] = dt;
        f[(1, 4)] = dt;
        f[(2, 5)] = dt;

        let q = Matrix6::identity() * 0.1; 

        self.state = f * self.state;
        self.covariance = f * self.covariance * f.transpose() + q;
    }

    fn update(&mut self, pulsars: &[Pulsar], measured_delays: &[f64], variances: &[f64]) -> Result<()> {
        let n = pulsars.len();
        
        let z = DVector::from_iterator(n, measured_delays.iter().map(|&d| -d * C_LIGHT));

        let mut h_rows = Vec::with_capacity(n * 6);
        for p in pulsars {
            h_rows.extend_from_slice(&[p.dir.x, p.dir.y, p.dir.z, 0.0, 0.0, 0.0]);
        }
        let h = DMatrix::from_row_slice(n, 6, &h_rows);

        let projected_measurement = &h * self.state;
        let y = &z - projected_measurement;

        let r = DMatrix::from_diagonal(&DVector::from_row_slice(variances));

        let s = &h * &self.covariance * h.transpose() + r;
        
        let s_inv = s.clone().try_inverse().ok_or(anyhow!("Singularidade na inversão"))?;
        let k = &self.covariance * h.transpose() * s_inv;

        self.state = &self.state + &k * y;

        let i = Matrix6::identity();
        self.covariance = (i - &k * &h) * self.covariance;

        Ok(())
    }
}

fn main() -> Result<()> {
    println!("--- XNAV: Inicializando Filtro de Kalman (OpenXNAV Core) ---");

    let mut pulsars = Vec::new();
    if let Ok(p) = Pulsar::from_par_file("J1937+21.par") { pulsars.push(p); }
    if let Ok(p) = Pulsar::from_par_file("J0437-4715.par") { pulsars.push(p); }
    
    if pulsars.is_empty() {
        println!("Aviso: Arquivos .par não encontrados. Usando catálogo interno.");
        pulsars = vec![
            Pulsar::new("PSR B1937+21", 20.0, 30.0, 0.00155, 500.0),
            Pulsar::new("PSR J0437-4715", 70.0, -47.0, 0.00575, 800.0), 
            Pulsar::new("PSR J1824-2452", 276.0, -24.0, 0.00305, 300.0),
            Pulsar::new("PSR J2124-3358", 321.0, -33.0, 0.00493, 250.0),
        ];
    } else {
        println!("Sucesso: Carregados {} pulsares via arquivos .par", pulsars.len());
    }

    let mut true_pos = Vector3::new(149_600_000.0, 0.0, 0.0); 
    let mut true_vel = Vector3::new(0.0, 30.0, 0.0);

    let initial_guess = Vector3::new(149_600_100.0, -50.0, 50.0);
    let mut nav_system = KalmanFilter::new(initial_guess);

    let mut rng = rand::thread_rng();

    let mut file = File::create("trajectory.csv")?;
    writeln!(file, "t,true_x,true_y,true_z,est_x,est_y,est_z,error_pos,uncertainty")?;

    println!("{:<5} | {:<15} | {:<15} | {:<15}", "Step", "Erro (km)", "Sigma (km)", "Evento");
    println!("{}", "-".repeat(65));

    for t in 0..100 {
        let dt = 1.0;

        let event = if t == 40 {
            true_vel += Vector3::new(0.0, 2.0, 0.0);
            "MANOBRA [BURN]"
        } else {
            ""
        };

        true_pos += true_vel * dt;

        let mut measurements = Vec::new();
        let mut noise_variances = Vec::new();

        for p in &pulsars {
            let roemer_delay = -(true_pos.dot(&p.dir)) / C_LIGHT;
            
            let shapiro_delay = p.shapiro_delay(&true_pos);

            let perfect_delay = roemer_delay + shapiro_delay;

            let photon_times = p.simulate_photons(dt, &mut rng);
            
            let measured_toa_noise = if !photon_times.is_empty() {
                let pulse_std_dev = p.width * p.period; 
                let n = photon_times.len() as f64;
                let standard_error = pulse_std_dev / n.sqrt();
                
                Normal::new(0.0, standard_error).unwrap().sample(&mut rng)
            } else {
                0.0
            };

            measurements.push(perfect_delay + measured_toa_noise);
            
            let sigma = if !photon_times.is_empty() {
                 (p.width * p.period) / (photon_times.len() as f64).sqrt()
            } else {
                1.0e9
            };
            
            noise_variances.push((sigma * C_LIGHT).powi(2));
        }

        nav_system.predict(dt); 
        // Atualizei a assinatura do update para aceitar variâncias dinâmicas
        nav_system.update(&pulsars, &measurements, &noise_variances)?;

        let estimated_pos = nav_system.state.fixed_rows::<3>(0);
        let error = (estimated_pos - true_pos).norm();
        let uncertainty = nav_system.covariance.diagonal().sum().sqrt();

        if t % 5 == 0 || t == 40 {
            println!("{:<5} | {:<15.4} | {:<15.4} | {}", 
                t, error, uncertainty, event
            );
        }

        writeln!(file, "{},{},{},{},{},{},{},{},{}", 
            t, 
            true_pos.x, true_pos.y, true_pos.z,
            estimated_pos[0], estimated_pos[1], estimated_pos[2],
            error, uncertainty
        )?;
    }

    println!("\nSimulação concluída. dados exportados para 'trajectory.csv'.");
    println!("use Matplotlib ou Excel para visualizar os resultados.");

    Ok(())
}

