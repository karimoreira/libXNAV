use nalgebra::{Vector3, Unit};
use std::f64::consts::PI;
use anyhow::{Result, Context};

const C_LIGHT: f64 = 299_792.458; 

#[derive(Debug, Clone)]
struct Pulsar {
    id: String,
    ra: f64,
    dec: f64,
  
    period: f64,
}

impl Pulsar {

    fn new(id: &str, ra_deg: f64, dec_deg: f64, period_s: f64) -> Result<Self> {
        Ok(Self {
            id: id.to_string(),
            ra: ra_deg.to_radians(),
            dec: dec_deg.to_radians(),
            period: period_s,
        })
    }


    fn direction_vector(&self) -> Unit<Vector3<f64>> {
        let x = self.dec.cos() * self.ra.cos();
        let y = self.dec.cos() * self.ra.sin();
        let z = self.dec.sin();
        Unit::new_normalize(Vector3::new(x, y, z))
    }
}


#[derive(Debug)]
struct Spacecraft { 
    position: Vector3<f64>,
    velocity: Vector3<f64>,
}

impl Spacecraft {
    fn new(x: f64, y: f64, z: f64) -> Self {
        Self {
            position: Vector3::new(x, y, z),
            velocity: Vector3::zeros(), 
        }
    }
}

struct XnavEngine;

impl XnavEngine {
    fn calculate_pulse_delay(spacecraft: &Spacecraft, pulsar: &Pulsar) -> f64 {
        let n = pulsar.direction_vector();
        let r = spacecraft.position;
        let projection = r.dot(&n);
        -projection / C_LIGHT
    }
}

fn main() -> Result<()> {
    println!("--- Iniciando simulação XNAV (rust) ---");

    let psr_1937 = Pulsar::new("PSR B1937+21", 294.5, 21.5, 0.00155)?;

    println!("alvo de navegação: {} (direção vetorial calculada)", psr_1937.id);

    let nave = Spacecraft::new(42_000.0, 15_000.0, -5_000.0);

    println!("posição da nave (km): {:?}", nave.position);
    let delay = XnavEngine::calculate_pulse_delay(&nave, &psr_1937);
    println!("\n--- resultados da Telemetria ---");
    println!("vetor direção do pulsar: {:.4?}", psr_1937.direction_vector().into_inner());
    println!("atraso de tempo calculado (geometric delay): {:.9} segundos", delay);
    println!("isso significa que o pulso chega {:.2} ns {} do que no baricentro.", 
             delay.abs() * 1e9,
             if delay < 0.0 { "antes" } else { "depois" }
    );

    Ok(())
}

