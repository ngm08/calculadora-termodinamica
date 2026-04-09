# Calculadora de Propriedades Termodinâmicas (Peng-Robinson EOS)

Esta é uma aplicação web desenvolvida para calcular propriedades termodinâmicas de substâncias puras, utilizando a equação de estado de Peng-Robinson. O objetivo é fornecer uma ferramenta simples e rápida para análises de equilíbrio líquido-vapor (ELV).

---

##  Fundamentação Teórica

A equação de **Peng-Robinson (1976)** é amplamente utilizada na engenharia química e de petróleo devido à sua precisão na estimativa de densidades de líquidos e propriedades próximas ao ponto crítico.

A expressão geral da pressão é dada por:

$$P = \frac{RT}{v - b} - \frac{a(T)}{v(v + b) + b(v - b)}$$

Onde:
* **$P$**: Pressão
* **$R$**: Constante universal dos gases
* **$T$**: Temperatura absoluta
* **$v$**: Volume molar
* **$a(T)$**: Parâmetro de atração (dependente da temperatura)
* **$b$**: Parâmetro de co-volume (volume ocupado pelas moléculas)

### Parâmetros de Atratividade e Co-volume:

Os parâmetros são calculados a partir das propriedades críticas ($P_c$, $T_c$) e do fator acêntrico ($\omega$), que são valores pré-determinados e tabelados para cada substância:

$$a(T) = \alpha \frac{0.45724 R^2 T_c^2}{P_c}$$
$$b = \frac{0.07780 R T_c}{P_c}$$

---

Para acessar a calculadora, basta clicar no link a seguir: https://calculadora-termodinamica.streamlit.app/
