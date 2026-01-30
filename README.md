# OpenXNAX Sistema de navegação baseado em Pulsares - XNAV

![Status](https://img.shields.io/badge/Status-Experimental-orange)
![Language](https://img.shields.io/badge/Language-Rust-red)
![License](https://img.shields.io/badge/License-MIT-blue)

Navegação autônoma baseada em raios-x

O OpenXNAV é um framework open-source para simulação de navegação de espaçonaves utilizando Pulsares de Milissegundos (MSPs). O projeto implementa modelos de propagação orbital e filtragem estocástica, focando na precisão dos modelos de tempo relativísticos.

O objetivo é fornecer uma base verificável para o estudo de algoritmos de determinação de órbita via XNAV, comparável a implementações de referência na literatura aeroespacial, visando democratizar a pesquisa acadêmica em XNAV

A navegação por pulsares baseia-se na estabilidade rotacional de estrelas de nêutrons. O sistema opera comparando o Tempo de Chegada (TOA) observado dos pulsos de raios-X com modelos preditivos calculados no Baricentro do Sistema Solar (SSB).

Metodologia de Estimação:

1.  **Aquisição de Sinal:** Simulação estocástica de contagem de fótons baseada no fluxo do pulsar.
2.  **Correção Baricêntrica:** Transformação de coordenadas temporais para o *Barycentric Celestial Reference System* (BCRS).
3.  **Resíduos de Tempo:** Cálculo da diferença entre observação geométrica e efemérides.
4.  **Filtragem Recursiva:** Aplicação de Filtro de Kalman Estendido (EKF) para fusão sensorial e propagação de covariância.

Arquitetura do Simulador

*   **Geração de Cenário:** Constelação configurável de pulsares com parâmetros astrométricos (Ascensão Reta, Declinação, Período, Fluxo).
*   **Propagador Orbital:** Integração numérica de estado (Posição/Velocidade/Aceleração).
*   **Sensor Virtual:** Modelagem de incerteza de medição baseada no limite de Cramer-Rao.


Instalação e Execução

Rust (Cargo) toolchain.

Simulação:

 - cargo run


A execução gera um log de telemetria (trajectory.csv) contendo a evolução temporal dos estados reais vs. estimados e a matriz de covariância.



Atraso de Shapiro (Relatividade Geral)
Implementação do atraso temporal gravitacional causado pela curvatura do espaço-tempo nas proximidades do Sol:
$$ \Delta t_{Shapiro} \approx - \frac{2GM}{c^3} \ln(1 - \cos \theta) $$
Este termo é mandatório para navegação precisa além de 1 UA.

Incerteza Radiométrica (Cramer-Rao)
A matriz de ruído de medição ($R$) não é estática. A variância ($\sigma^2$) é derivada dinamicamente das propriedades físicas do pulsar e do tempo de integração ($T_{obs}$):
$$ \sigma_{TOA} \propto \frac{W}{F \cdot \sqrt{T_{obs}}} $$
Isso permite análises de sensibilidade sobre a área do detector e geometria da constelação.

Motor de Estimativa (EKF 6-DoF)
Correções Relativísticas (Atraso de Shapiro)
Modelo de Ruído dependente de Fluxo (Physics-based)
Leitor de Efemérides (.par files)
Simulação de Processo de Poisson para fótons individuais


Contribuições focadas em otimização numérica, novos propagadores orbitais ou refinamento dos modelos de tempo são bem vindas via pull requests
