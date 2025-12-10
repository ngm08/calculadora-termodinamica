import streamlit as st
import math

# --- CONFIGURAÇÃO DA PÁGINA ---
# Mudei para "wide" para facilitar o alinhamento lado a lado, 
# mas funciona em "centered" também se preferir.
st.set_page_config(page_title="Calculadora Termodinâmica", layout="wide") 

st.title("⚗️ Calculadora Termodinâmica")
st.markdown("Cálculo de propriedades termodinâmicas com base na equação de estado de **Peng-Robinson**.")

# --- DADOS E CONSTANTES ---
R = 8.314e-5 #bar.m³/K.mol
Tref= 273.15 
Pref= 1.01325
dados = {
            'Acetileno': {
                'pc': 61.39, 'tc': 308.3, 'w': 0.187,
                'const_Cp': [2.410,10.926e-3,-0.255e-5,-0.790e-8,0.524e-11]
            },
            'Água': {
                'pc': 221.2, 'tc': 647.3, 'w': 0.344,
                'const_Cp': [4.395,-4.186e-3,1.405e-5,-1.564e-8,0.632e-11]
            },
            'Benzeno': {
                'pc': 48.98, 'tc': 562.2, 'w': 0.211,
                'const_Cp': [3.551,-6.184e-3,14.365e-5,-19.807e-8,8.234e-11]
            },
            'Ciclo-hexano': {
                'pc': 40.805, 'tc': 553.600, 'w': 0.20926,
                'const_Cp': [3.874, -0.909e-3, 14.902e-5, -19.907e-8, 8.011e-11]
            },
            'Dióxido de Carbono': {
                'pc': 73.82, 'tc': 304.2, 'w': 0.228,
                'const_Cp': [3.259,1.356e-3,1.502e-5,-2.374e-8,1.056e-11]
            },
            'Etano': {
                'pc': 48.8, 'tc': 305.4, 'w': 0.099,
                'const_Cp': [4.178,-4.427e-3,5.660e-5,-6.651e-8,2.487e-11]
            },
            'Etileno': {
                'pc': 50.32, 'tc': 282.4, 'w': 0.085,
                'const_Cp': [4.221,-8.782e-3,5.795e-5,-6.729e-8,2.511e-11]
            },
            'Hidrogênio': {
                'pc': 12.97, 'tc': 33.3, 'w': -0.215,
                'const_Cp': [2.883,3.681e-3,-0.772e-5,0.692e-8,-0.213e-11]
            },
            'Metano': {
                'pc': 46.04, 'tc': 190.6, 'w': 0.011,
                'const_Cp': [4.568,-8.975e-3,3.631e-5,-3.407e-8,1.091e-11]
            },
            'n-Butano': {
                'pc': 37.97, 'tc': 425.2, 'w': 0.193,
                'const_Cp': [5.547,5.536e-3,8.057e-5,-10.571e-8,4.134e-11]
            },
            'n-Hexano': {
                'pc': 30.12, 'tc': 507.4, 'w': 0.305,
                'const_Cp': [8.831,-0.166e-3,14.302e-5,-18.314e-8,7.124e-11]
            },
            'n-Pentano': {
                'pc': 33.69, 'tc': 469.7, 'w': 0.249,
                'const_Cp': [7.554,-0.368e-3,11.846e-5,-14.939e-8,5.753e-11]
            },
            'Nitrogênio': {
                'pc': 33.94, 'tc': 126.1, 'w': 0.04,
                'const_Cp': [3.539,-0.261e-3,0.007e-5,0.157e-8,-0.099e-11]
            },
            'Oxigênio': {
                'pc': 50.43, 'tc': 154.6, 'w': 0.022,
                'const_Cp': [3.630,-1.794e-3,0.658e-5,-0.601e-8,0.179e-11]
            },
            'Propano': {
                'pc': 42.49, 'tc': 369.8, 'w': 0.152,
                'const_Cp': [3.847,5.831e-3,6.013e-5,-7.893e-8,3.079e-11]
            },
            'Propileno': {
                'pc': 46.13, 'tc': 364.8, 'w': 0.142,
                'const_Cp': [3.834,3.893e-3,4.688e-5,-6.013e-8,2.283e-11]
            },
            'Tolueno': {
                'pc': 41.09, 'tc': 591.8, 'w': 0.264,
                'const_Cp': [3.866,3.558e-3,13.356e-5,-18.659e-8,7.690e-11]
            }
        }

# --- FUNÇÕES DE CÁLCULO ---
def alfa(tr, w):
    k = 0.37464 + 1.54226*w-0.26992*w**(2)
    alfa = (1 + k*(1 - (tr)**0.5))**2
    return alfa,k

def AB(p, t, pc, tc, w):
    tr = t/tc
    alpha,k = alfa(tr, w)
    a = 0.45724*(R**(2)*tc**(2)/pc)*alpha
    b = 0.07780*(R*tc/pc)
    A = (a*p)/(R*t)**(2)
    B = (b*p)/(R*t)
    return A, B, a, b, alpha, k

def f(Z, A, B):
    return Z**(3)-(1-B)*Z**(2)+(A - 3*B**2 - 2*B)*Z - (A*B - B**2 - B**3)

def df(Z, A, B):
    return 3*Z**(2)-2*(1-B)*Z+A-3*B**2 - 2*B

def nrz(f, df, Z_chute, A, B, max_iter=200):
    Z_old = Z_chute
    erro = 1.0
    tol = 1e-8
    i = 0
    while i < max_iter:
        den = df(Z_old,A,B)
        if den == 0: return None 
        Z_new = Z_old - f(Z_old,A,B)/den
        erro = abs(Z_new - Z_old)
        Z_old = Z_new 
        if erro < tol:
            return Z_old
        i += 1
    return None

def Cp_T(T, Tref, Cp_coefs):
    return R*(Cp_coefs[0]*math.log(T/Tref) + Cp_coefs[1]*(T - Tref) + (Cp_coefs[2]/2)*(T**2 - Tref**2) + (Cp_coefs[3]/3)*(T**3 - Tref**3) + (Cp_coefs[4]/4)*(T**4 - Tref**4))

def h_res(Z, T, P, Pc, Tc, w):
    A, B, a, b, alpha, k = AB(P, T, Pc, Tc, w)
    da = -0.45724*((R**2)*(Tc**2)/Pc)*k*((alpha/(T*Tc))**(1/2))
    arg = (Z+(1+(2**(1/2)))*B)/(Z+(1-(2**(1/2)))*B)
    if arg <= 0: arg = 1e-10 
    hr = (R*T*(Z-1)+((T*da-a)/(2*b*(2**(1/2))))*math.log(arg))*1e5
    return hr

def delta_h(Tref, T, Cp_coefs):
    ans = (R*(Cp_coefs[0]*(T-Tref)+Cp_coefs[1]*((T**2-Tref**2))/2+Cp_coefs[2]*((T**3-Tref**3))/3+Cp_coefs[3]*((T**4-Tref**4))/4+Cp_coefs[4]*((T**5-Tref**5))/5))*1e5
    return ans

def h_calc(Z, T, P, Tref, Pref, Pc, Tc, Cp_coefs, w):
    ans = h_res(Z, T, P, Pc, Tc, w) + delta_h(Tref, T, Cp_coefs) - h_res(Z, Tref, P, Pc, Tc, w)
    return ans

def s_res(Z, T, P, Pc, Tc, w):
    A, B, a, b, alpha, k = AB(P, T, Pc, Tc, w)
    da= -0.45724*((R**2)*(Tc**2)/Pc)*k*((alpha/(T*Tc))**(1/2))
    arg_log1 = Z-B
    if arg_log1 <= 0: arg_log1 = 1e-10
    arg_log2 = (Z+(1+(2**(1/2)))*B)/(Z+(1-(2**(1/2)))*B)
    if arg_log2 <= 0: arg_log2 = 1e-10
    
    sr = (R*math.log(arg_log1)+((da)/(2*b*(2**(1/2))))*math.log(arg_log2))*1e5
    return sr

def delta_s(Tref, T, P, Pref, Cp_coefs):
    ans = ((Cp_T(T, Tref, Cp_coefs)) - R*math.log(P/Pref))*1e5
    return ans

def s_calc(Z, T, P, Tref, Pref, Pc, Tc, Cp_coefs, w):
    ans = s_res(Z, T, P, Pc, Tc, w) + delta_s(Tref, T, P, Pref, Cp_coefs) - s_res(Z, Tref, Pref, Pc, Tc, w)
    return ans

def calc_Z(A, B, P, Pc):
    raizes = []
    for chute in [i/20 for i in range(1, 21)]: 
        raiz = nrz(f, df, chute, A, B)
        if raiz == None:
            continue
        if not any(abs(raiz - r) < 1e-8 for r in raizes):
            if raiz > 0: 
                raizes.append(raiz)

    Z_lista = sorted(raizes)
    
    if not Z_lista: return 1.0 

    if len(Z_lista) == 1:
        return Z_lista[0]
    else:
        return Z_lista

def turbina(P_in, T_in, P_out, eff, taxa, Pref, Tref, Composto, tol=1e-5):
    Pc, Tc, w, Cp_coefs = dados[Composto].values()
    A, B, *_ = AB(P_in, T_in, Pc, Tc, w)
    Z_in = calc_Z(A, B, P_in, Pc)
    
    if isinstance(Z_in, list):
        if P_in > Pc: Z_in = Z_in[0] 
        else: Z_in = Z_in[-1] 

    entalpia_in = h_calc(Z_in, T_in, P_in, Tref, Pref, Pc, Tc, Cp_coefs, w)
    entropia_in = s_calc(Z_in, T_in, P_in, Tref, Pref, Pc, Tc, Cp_coefs, w)
    
    delta_T = 1e-4
    T_out_isentropica = T_in
    for _ in range(500): 
        A_out, B_out, *_ = AB(P_out, T_out_isentropica, Pc, Tc, w)
        Z_out_s = calc_Z(A_out, B_out, P_out, Pc)
        if isinstance(Z_out_s, list): Z_out_s = Z_out_s[-1] 

        entropia_out = s_calc(Z_out_s, T_out_isentropica, P_out, Tref, Pref, Pc, Tc, Cp_coefs, w)
        
        T_corrigida = T_out_isentropica + delta_T
        A_plus, B_plus, *_ = AB(P_out, T_corrigida, Pc, Tc, w)
        Z_plus = calc_Z(A_plus, B_plus, P_out, Pc)
        if isinstance(Z_plus, list): Z_plus = Z_plus[-1]
        
        entropia_plus = s_calc(Z_plus, T_corrigida, P_out, Tref, Pref, Pc, Tc, Cp_coefs, w)
        
        diff = entropia_out - entropia_in
        if abs(diff) < tol: break
        
        df = (entropia_plus - entropia_out) / delta_T
        if df == 0: break
        dT = diff/df
        T_out_isentropica -= dT
    
    entalpia_isentropica = h_calc(Z_out_s, T_out_isentropica, P_out, Tref, Pref, Pc, Tc, Cp_coefs, w)
    
    potencia = -(eff * (entalpia_isentropica - entalpia_in)*taxa)/1000
    entalpia_out_real = (eff * (entalpia_isentropica - entalpia_in)) + entalpia_in
    
    T_out = T_out_isentropica
    for _ in range(500):
        A_real, B_real, *_ = AB(P_out, T_out, Pc, Tc, w)
        Z_real = calc_Z(A_real, B_real, P_out, Pc)
        if isinstance(Z_real, list): Z_real = Z_real[-1]

        entalpia_out = h_calc(Z_real, T_out, P_out, Tref, Pref, Pc, Tc, Cp_coefs, w)
        
        diff_h = entalpia_out - entalpia_out_real
        if abs(diff_h) < tol: break
        
        T_plus = T_out + delta_T
        A_plus, B_plus, *_ = AB(P_out, T_plus, Pc, Tc, w)
        Z_plus = calc_Z(A_plus, B_plus, P_out, Pc)
        if isinstance(Z_plus, list): Z_plus = Z_plus[-1]
        
        entalpia_plus = h_calc(Z_plus, T_plus, P_out, Tref, Pref, Pc, Tc, Cp_coefs, w)
        
        df = (entalpia_plus - entalpia_out) / delta_T
        if df == 0: break
        dT = diff_h/df
        T_out -= dT

    entropia_out_real = s_calc(Z_real, T_out, P_out, Tref, Pref, Pc, Tc, Cp_coefs, w)
    entropia_produzida = (entropia_out_real - entropia_in)*taxa
    
    return potencia, T_out, entropia_produzida

# --- INTERFACE PRINCIPAL ---

st.header("1. Propriedades de Estado")
st.markdown("Cálculo de volume, entalpia e entropia para uma condição específica.")

# Distribuição das colunas da seção 1
col_prop1, col_prop2, col_prop3, col_prop4 = st.columns(4)

with col_prop1:
    composto_prop = st.selectbox("Substância:", list(dados.keys()), index=0)
with col_prop2:
    P_prop = st.number_input("Pressão (bar):", value=10.0, min_value=0.1, step=1.0, key="p_prop")
with col_prop3:
    T_prop = st.number_input("Temperatura (K):", value=350.0, min_value=100.0, step=1.0, key="t_prop")
with col_prop4:
    st.write("") # Espaçamento para alinhar botão
    st.write("")
    calc_btn = st.button("Calcular", use_container_width=True)

if calc_btn:
    Pc, Tc, w, Cp_coefs = dados[composto_prop].values()
    A, B, *_ = AB(P_prop, T_prop, Pc, Tc, w)
    Z = calc_Z(A, B, P_prop, Pc)
    
    if isinstance(Z, list):
        st.info(f"Mais de uma fase detectada. Valores possíveis para **Z: {[round(r, 4) for r in Z]}**")
        # Simplificação: se P > Pc assume denso, senão vapor (lógica comum em PR simples)
        if P_prop > Pc:
            Z_final = Z[0]
            st.write("Fase considerada: **Líquido Comprimido** (Considera o menor valor de Z)")
        else:
            Z_final = Z[-1]
            st.write("Fase considerada: **Vapor** (Considera o maior valor de Z)")
    else:
        Z_final = Z
        st.write(f"Fase única detectada. **Z = {Z_final:.4f}**")
    
    vm = (Z_final * R * T_prop / P_prop)*1e3
    h_val = h_calc(Z_final, T_prop, P_prop, Tref, Pref, Pc, Tc, Cp_coefs, w)
    s_val = s_calc(Z_final, T_prop, P_prop, Tref, Pref, Pc, Tc, Cp_coefs, w)

    res_c1, res_c2, res_c3 = st.columns(3)
    res_c1.success(f"**Volume Molar:**\n\n {vm:.4f} L/mol")
    res_c2.success(f"**Entalpia:**\n\n {h_val:.2f} J/mol")
    res_c3.success(f"**Entropia:**\n\n {s_val:.4f} J/mol.K")

st.markdown("---") 

# --- PARTE 2: SIMULAÇÃO DE TURBINA ---
st.header("2. Simulação de Turbina")
st.markdown("Insira as informações de entrada e saída de sua turbina.")

# ALINHAMENTO FEITO AQUI: 3 Inputs na esquerda, 3 Inputs na direita
t_col1, t_col2 = st.columns(2)

with t_col1:
    st.subheader("Entrada")
    composto_turb = st.selectbox("Substância:", list(dados.keys()), key="tb_comp")
    P_in_t = st.number_input("Pressão de Entrada (bar)", value=60.0, step=1.0)
    T_in_t = st.number_input("Temperatura de Entrada (K)", value=650.0, step=1.0)

with t_col2:
    st.subheader("Saída e Parâmetros")
    taxa = st.number_input("Taxa de fluxo (mol/s)", value=20.0, step=1.0) # Movido para cá para equilibrar
    P_out_t = st.number_input("Pressão de Saída (bar)", value=10.0, step=1.0)
    eff = st.slider("Eficiência Isentrópica", 0.0, 1.0, 0.8)

if st.button("Simular Turbina", type="primary", use_container_width=True):
    with st.spinner(f'Calculando expansão de {composto_turb}...'):
        try:
            pot, T_out_res, s_ger = turbina(P_in_t, T_in_t, P_out_t, eff, taxa, Pref, Tref, composto_turb)
            
            st.divider()
            r1, r2, r3 = st.columns(3)
            r1.metric("⚡ Potência", f"{pot:.2f} kW")
            r2.metric("🌡️ Temperatura de Saída", f"{T_out_res:.2f} K")
            r3.metric("🔥 Entropia Gerada", f"{s_ger:.4f} kJ/K")
            
        except Exception as e:
            st.error(f"Erro no cálculo: {e}. Verifique se as condições de entrada são fisicamente possíveis.")