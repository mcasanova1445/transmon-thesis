---
title:
- Diseño y simulación de procesadores cuánticos que implementen algoritmos cuánticos de busqueda
author:
- Miguel Casanova
header-includes:
   - \usepackage{physics}
   - \usepackage{mathtools}
   - \usepackage{physics}
   - \usepackage{qcircuit}
   - \usepackage{graphicx}
   - \usepackage{float}
   - \usepackage{dsfont}
   - \usepackage{tkz-graph}
   - \usetikzlibrary{arrows}
   - \usepackage{ragged2e}
theme:
- Copenhagen
---


# Estructura de la presentación

# Objetivos

# Información cuántica

# Espacios de Hilbert

Un espacio de Hilbert es un espacio lineal real o complejo con un producto interno que también define un espacio normado completo.

- Producto interno
- Normado
- Completo: Todas las sucesiones de Cauchy convergen fuertemente.
    * Sucesiones de Cauchy: Aquellas en las que el ordenamiento no afecta la convergencia.
- Separable: Tiene bases contables.

# Operadores

- Operadores hermíticos $U = U^\dagger$
    * Autovalores reales
    * Diagonal real
    * Diagonalizable

- Operadores unitarios: $U U^\dagger = \mathds{1}$
    * Determinante de módulo igual a la unidad
    * Preserva normas y trazas
    * Diagonalizable

# Estados cuánticos

- Estados puros: $\ket{\psi}$
- Estados mixtos: $\rho$

# Sistemas multipartitos

$\mathcal{H} = \mathcal{H}_1 \mathcal{H}_2$

$\rho_B = Tr_A(\rho_{AB})$

# Qubits y esfera de Bloch

Los qubits son la unidad básica de información cuántica.

$\ket{0} = \begin{pmatrix}1 \\ 0\end{pmatrix}, \quad \ket{1} = \begin{pmatrix}0 \\ 1\end{pmatrix}$

$\ket{\psi} = \alpha \ket{0} + \beta \ket{1}, \quad \alpha, \beta \in \mathds{C}$

Esfera de Bloch:

# Ecuación de Schrödinger

Esta ecuación describe la evolución de los estados cuánticos puros.

$i \hbar \frac{d}{dt} \ket{\psi(t)} = \hat{H} \ket{\psi(t)}$

Para estados mixtos se utiliza la ecuación de Liouville-von Neumann.

$i \hbar \dot{\rho}(t) = [\hat{H}, \rho(t)]$

# Postulados de la mecánica cuántica

# Matrices de Pauli

$\sigma_x = \begin{pmatrix}0 & 1 \\ 1 & 0\end{pmatrix}$

$\sigma_y = \begin{pmatrix}0 & -i \\ i & 0\end{pmatrix}$

$\sigma_z = \begin{pmatrix}1 & 0 \\ 0 & -1\end{pmatrix}$

# Circuitos cuánticos

Circuito de ejemplo

# Compuertas cuánticas

Compuertas con matrices

# Criterios de Di Vincenzo

- Un sistema físico escalable con qubits caracterizados.
- La habilidad de inicializar el estado de los qubits en un estado fiducial simples.
- Tiempos de coherencia relevantes largos.
- Un conjunto universal de compuertas cuánticas.
- La capacidad de medir qubits en específico.

# Fidelidad

La fidelidad es una medida de distancia entre estados cuánticos. Dos estados idénticos tendrán una fidelidad igual a 1, mientras que dos estados ortogonales tendrán una fidelidad de 0.

$F(\rho, \sigma) = Tr(\sqrt{\rho^{1/2} \sigma \rho^{1/2}})$

# Medidas proyectivas

# Sistemas cuánticos abiertos

En la mecánica cuántica de sistemas abiertos con evolución markoviana, la ecuación de Schrödinger toma la siguiente forma más general, conocida como Lindbladiano.

$\dot{\rho}(t) = -i [\hat{H}, \rho(t)] + \sum_k \gamma_k [V_k \rho(t) V_k^\dagger - \frac{1}{2} \{V_k^\dagger V_k, \rho(t)\}]$

# Superconductividad

# Teorías BCS

# Efecto Josephson

\justify
Se tienen dos placas superconductoras A y B, separadas por un aislante.
Las funciones de onda de las placas superconductoras son:
$\psi_A = \sqrt{\rho_1} e^{i \phi_1}, \psi_B = \sqrt{\rho_2} e^{i \phi_2}$ 

\justify
En principio, no debería haber conducción entre ambas placas.
Sin embargo, ese no es el caso. Por el efecto tunel, una supercorriente (corriente
sin disipación) de pares de Cooper (pares de electrones con spines opuestos)
pueden pasar de una placa a la otra sin disipación.

$$V_J = \frac{\hbar}{2e} \frac{d\delta}{dt}$$

$$I_J = I_0 \sin(\delta)$$

\justify
Donde $\delta=\phi_2-\phi_1$ es la diferencia de fase entre las dos placas superconductoras.

# Efecto Josephson DC y AC

* Efecto Josephson DC

\justify
Si las placas se encuentran sin alimentación, entonces correrá una
supercorriente constante a través de ellas.

* Efecto Josephson AC

\justify
Si las placas se alimentan con un voltaje DC externo, entonces la diferencia
de fase entre ellas variará linealmente con el tiempo y habrá una corriente AC
a través de ellas.

# Energía e inductancia de Josephson

$$E_J = \int I_0 \sin(\delta) \frac{\hbar}{2e} \frac{d\delta}{dt} dt
 = \frac{\hbar I_0}{2e} \int \sin(\delta) dt
 = - \frac{\hbar I_0}{2e} \cos(\delta)$$

$\frac{dI_J}{dt} = I_0 \cos(\delta) \frac{d \delta}{dt} = I_0 \cos(\delta)
 \frac{2e}{\hbar} V_J$

$$L_J = \frac{\hbar}{2e I_0 \cos(\delta)}$$

$E_C = \frac{(2e)^2}{2C}$

$E_L = \frac{\hbar^2}{(2e)^2L}$

# Circuito LCJ

\justify
Primero, reescribimos el Hamiltoniano del circuito LC en términos de la
cantidad de pares de Cooper y de la diferencia de fase en los extremos
del inductor, en lugar de la carga y el flujo

$$\hat{H} = \frac{(2e)^2}{2C} \hat{n}^2 + \frac{\hbar^2}{(2e)^2L} 
 \frac{\hat{\delta^2}}{2}$$

De aquí $\hat{q}=2e\hat{n}$ y $\hat{\phi}=\frac{\hbar}{2e}\hat{\delta}$

Ahora introducimos el término de la unión Josephson

$\hat{H} = \frac{(2e)^2}{2C} (\hat{n}-n_g)^2 + \frac{\hbar^2}{(2e)^2L} 
 \frac{\hat{(\delta-\delta_e)^2}}{2} - \frac{\hbar I_0}{2e} \cos(\delta)
 = E_C (\hat{n}-n_g)^2 + E_L \frac{(\hat{\delta}-\delta_e)^2}{2}
 - E_{J0} \cos(\hat{\delta})$

# Ecuación de Schrödinger del circuito LCJ

\justify
Para describir el sistema en términos de la ecuación de Schrödinger en
función de la diferencia de fase $\delta$, se introduce $\hat{n}=-i\hbar
\frac{\partial}{\partial \phi}$

$$E_C (-i\hbar \frac{\partial}{\partial\phi}-n_g)^2+U(\phi))\psi = E \psi$$
$$U(\phi) = -E_{J0} \cos(\phi)+E_L \frac{(\phi-\phi_e)^2}{2}$$

# Arquetipos de qubits superconductores

\justify
* Qubit de carga: Si $E_L$ tiende a cero, la carga almacenada en la isla 
    superconductora entre el capacitor y  la unión Josephson se puede usar 
    como qubit. El potencial de este tipo de qubit es de forma de coseno.

\justify
* Qubit de flujo: Si $E_L$ es comparable con $E_{J0}$, el flujo a través 
    del lazo formado por el inductor y la unión Josephson se puede usar como 
    qubit. El potencial de este tipo de qubit es de forma cuártica.

\justify
* Qubit de fase: Si se polariza la unión Josephson con una fuente de 
    corriente, la fase en ambos extermos de la unión Josephson se puede 
    usar como qubit. El potencial de este tipo de qubit es de forma cúbica.

# Qubits de carga

$[\hat{\delta},\hat{n}]=i \implies e^{\pm i \hat{\delta}} \ket{n} = \ket{n\pm 1}$

$\hat{H} = E_C (\hat{n}-n_g)^2 - E_{J0} \cos(\hat{\delta}) =
 E_C (\hat{n}-n_g)^2 - E_{J0} (e^{i \hat{\delta}} + e^{-i \hat{\delta}}) =
 \sum(E_C(N-N_g)^2 \ketbra{n}{n} - \frac{E_{J0}}{2}(\ketbra{n}{n+1}+
 \ketbra{n+1}{n}))$

# Caja de pares de Cooper 

Dimensiones típicas de la isla: 1000nm x 50nm x 20nm

![Circuito de una caja de pares de Cooper](Avance1/cooperpairbox.png){#id .class width=45%} ![Niveles de energía](Avance1/cooperenergy.png){#id .class width=45%}

# Transmon

Intercambiamos anarmonicidad por independencia de $n_g$

![Circuito de un transmon](Avance1/transmon.png){#id .class width=45%} ![Niveles de energía](Avance1/transmonenergy.png){#id .class width=45%}

# Modelo de Jaynes-Cummings

$$\hat{H} = \hat{H}_c + \hat{H}_q + \hat{H}_{qc} = \hbar \omega_c (a^\dag a + 
  \frac{1}{2}) + \frac{1}{2} \hbar \omega_q \sigma_z + \hbar g \sigma_x (a+a^\dag)$$

\justify
De ahora en adelante $\hbar = 1$ y despreciaré los términos constantes, pues 
sólo contribuyen en fases globales a la evolución del sistema.

# Aproximación de onda rotacional

$$\hat{H}_{qc} = \hat{H}_{qc}^{JC} + \hat{H}_{qc}^{AJC} = g(a \sigma_+ + 
  a^\dag\sigma_-) + g(a^\dag \sigma_+ + a \sigma_-)$$

$$\hat{H} = \hat{H}_c + \hat{H}_q + \hat{H}_{qc} = \omega_c a^\dag a +
  \frac{1}{2} \omega_q \sigma_z + g(a \sigma_+ + a^\dag \sigma_-)$$

# Hamiltoniano multiquibit

$$\hat{H} = \hat{H}_q + \hat{H}_{qc} = \frac{1}{2} \sum\limits_i \omega_{qi} 
  \sigma_{zi} + \sum\limits_i g_i (a \sigma_{+ i} + a^\dagger \sigma_{- i})$$

# Pulsos de microondas

$$\hat{H}_d = \sum\limits_k (a+a^\dagger) (\xi_k e^{-i\omega_d^{(k)}t} + 
  \xi_k^*e^{i\omega_d^{(k)}t})$$

RWA: $$\hat{H}_d=\sum\limits_k a\xi_k^*e^{i\omega_d^{(k)}t}+
  a^\dagger\xi_ke^{-i\omega_d^{(k)}t}$$

# Régimen rotacional del pulso

\justify
Trabajando con un sólo modo a la vez, se aplica la siguiente transformación $U(t) 
= exp[-i \omega_d t(a^\dagger a + \sum\limits_i \sigma_{z i})]$ para entrar en 
el régimen rotacional del pulso de control.

$$\hat{H} = U^\dagger (\hat{H}_{syst} + \hat{H}_d) U - i U^\dagger \dot{U}$$
$$ \hat{H} = \Delta_c a^\dagger a + \frac{1}{2} \sum\limits_i \Delta_{qi} 
  \sigma_{zi} + \sum\limits_i g_i (a \sigma_{+ i} + a^\dagger \sigma_{- i}) + 
  (a\xi^*e^{i\omega_d t}+a^\dagger\xi e^{-i\omega_d t})$$

$\Delta_c = \omega_c - \omega_d \qquad \quad \Delta_{qi} = \omega_{qi} - \omega_d$

# Efecto del pulso sobre el qubit

\justify
Luego se aplica el operador de desplazamineto $D(\alpha) = exp[\alpha a^\dagger - 
\alpha^* a]$ sobre el campo $a$ con $\dot{\alpha} = -i \Delta_c \alpha -i \xi 
e^{-i \omega_d t}$ para eliminar el efecto directo del pulso sobre la cavidad.

$$\hat{H} = D^\dagger (\alpha) \hat{H}_{old} D(\alpha) -i D^\dagger(\alpha) \dot{D}(\alpha)$$

$$\hat{H} = \Delta_c a^\dagger a + \frac{1}{2} \sum\limits_i \Delta_{qi} 
  \sigma_{zi} + \sum\limits_i g_i (a \sigma_{+i} + a^\dagger \sigma_{-i})$$
$$ + \sum\limits_i g_i (\alpha \sigma_{+i} + \alpha^* \sigma_{-i}) - \Delta_c
  \alpha \alpha^* $$

El término $-\Delta_c \alpha \alpha^*$ se desprecia, ya que sólo representa 
una fase global en la evolución del sistema.

# Régimen dispersivo

\justify
Finalmente, aplicamos la transformación $U = exp[\sum\limits_i \frac{g_i}
{\Delta_i} (a^\dagger \sigma_{-i} - a \sigma_{+i})]$, donde $\Delta_i = 
\omega_{qi} - \omega_c$ y realizamos la aproximación de segundo grado sobre
los términos $\frac{g_i}{\Delta_i} \ll 1$.

$$\hat{H} = U^\dagger \hat{H}_{old} U$$
$$\hat{H} \approx \tilde{\Delta}_c a^\dagger a + \frac{1}{2} \sum\limits_i
  \tilde{\Delta}_{qi} \sigma_{zi} + \sum\limits_i (\Omega_i \sigma_{+i} +
  \Omega_i^* \sigma_{-i})$$
$$+ \sum\limits_{i \neq j} \frac{g_i g_j}{2 \Delta_i} 
  (\sigma_{-i} \sigma_{+j}+\sigma_{+i} \sigma_{-j})$$

$\tilde{\Delta}_c = (\omega_c + \sum\limits_i \chi_i \sigma_{zi}) - \omega_d
 \qquad
 \tilde{\Delta}_{qi} = (\omega_{qi} + \chi_i) - \omega_d
 \qquad
 \chi_i = \frac{g_i^2}{\Delta_i}$

# Rotaciones X-Y

\justify
Tomando $\Omega(t) = \Omega^x(t) \cos(\omega_d t) + \Omega^y \sin(\omega_d t)$,
donde $\omega_d$ es igual a la frecuencia de resonancia de uno de los qubits
logramos rotaciones sobre los ejes X e Y. Las amplitudes de estas rotaciones
vienen dadas por $\int_0^{t_0} \Omega^x(t) dt$ y $\int_0^{t_0} \Omega^y(t)
dt$, respectivamente, donde $t_0$ es la duración del pulso.

$$\hat{H} \approx \tilde{\Delta}_c a^\dagger a + \frac{1}{2} \tilde{\Delta}_q 
  \sigma_z + \frac{1}{2} (\Omega^x(t) \sigma_x + \Omega^y(t) \sigma_y)$$

# Compuerta de entrelazamiento

Ejemplo con sólo dos qubits

$$\hat{H} \approx \frac{1}{2} \tilde{\Delta}_{q_1} \sigma_{z_1} + 
  \frac{1}{2} \tilde{\Delta}_{q_2} \sigma_{z_2} + 
  \frac{g_1 g_2 (\Delta_1 + \Delta_2)}{2 \Delta_1 \Delta_2} 
  (\sigma_{-_1} \sigma_{+_2} + \sigma_{+_1} \sigma_{-_2})$$

Variando la frecuencia de resonacia de los qubit, se puede variar el 
acoplamiento entre estos. 

# Simulador

# Parámetros de los sistemas simulados

# Compuertas compuestas

# Algoritmo de Grover

\begin{align}
    \Pi_\mathcal{W} &= \sum_k \ketbra{\omega_k} \\
    \ket{\psi} &= \sin(\theta) \ket{\psi_1} + \cos(\theta) \ket{\psi_0}
\end{align}

\begin{align}
    U_\mathcal{W} &= 2 \mathds{1} - 2 \Pi_\mathcal{W}\\
    U_\psi &= 2 \ketbra{\psi} - \mathds{1}
\end{align}


# Variaciones del algoritmo de Grover

# Simulaciones del algoritmo de Grover

# Transformada cuántica de Fourier

\begin{equation}
    QFT^\dagger = \frac{1}{\sqrt{N}} \sum_x \sum_k e^{-2 \pi i k x / N} \ketbra{x}{k}
\end{equation}

# Estimación de fase

# Estimación de orden

# Algoritmo de Shor

# Simulaciones del algoritmo de Shor

# PageRank

# Caminatas cuánticas de Szegedy

# PageRank cuántico

# Circuitos de Loke

# Simulaciones del algoritmo de PageRank cuántico

# Conclusiones






































































