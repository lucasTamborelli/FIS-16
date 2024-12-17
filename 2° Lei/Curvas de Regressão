import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import stats
import statistics

print('quantidade de experimentos:')  # 12
n = input()
n = int(n)

print('deslocamento 1 e 2 (m)')  # 0.336 e 0.686
ds1 = input()
ds1 = float(ds1)
ds2 = input()
ds2 = float(ds2)

f = open("dadosEXP3.txt", "r")
t_1 = 0.0
t_2 = 0.0
m = 0.0
tempos, massas = [], []
st1 = [0, 0, 0]
st2 = [0, 0, 0]
sigma_t = []

for linha, valor in enumerate(f.readlines()):
    valor = float(valor.strip())

    if linha % 7 == 0:

        if linha != 0:
            tempos.append(t_1)
            tempos.append(t_2)
            sigma_t.append(statistics.stdev(st1))
            sigma_t.append(statistics.stdev(st2))
            massas.append(m)
            massas.append(m)
            t_1 = 0
            t_2 = 0

        m = valor / 1000

    elif linha % 7 == 1:
        t_1 = t_1 + valor / 3
        st1[0] = valor

    elif linha % 7 == 2:
        t_2 = t_2 + valor / 3
        st2[0] = valor

    elif linha % 7 == 3:
        t_1 = t_1 + valor / 3
        st1[1] = valor

    elif (linha % 7 == 4):
        t_2 = t_2 + valor / 3
        st2[1] = valor

    elif (linha % 7 == 5):
        t_1 = t_1 + valor / 3
        st1[2] = valor

    elif (linha % 7 == 6):
        t_2 = t_2 + valor / 3
        st2[2] = valor

temposQuadrado = np.array(list(map(lambda x: x * x, tempos)))
deslocamentos = np.array([ds1, ds2] * n)
y = 2 * np.divide(deslocamentos, temposQuadrado)

invM = np.array(1 / np.array(massas))  # vetor com o inverso das massas


# CALCULO DO ERRO FINAL DE Y:
# considerando x (inverso das Massas) insento de erro

sigma_t = np.array(sigma_t)
tempos = np.array(tempos)
sigma_y = np.array(np.sqrt((2 * 0.0005 / (tempos) ** 2) ** 2 + (4 * sigma_t * deslocamentos / (tempos) ** 3) ** 2))
# print(sigma_y)
sigma_invM = np.array(0.00001 * (invM) ** 2)


def modelo(x, a, b):
    return a * x + b


popti, pcovi = curve_fit(modelo, invM, y, sigma=sigma_y, absolute_sigma=True)
a_i, b_i = popti
print('a para o calculo do ErroFinal', a_i)

erroFinal = np.sqrt(sigma_y ** 2 + a_i * a_i * (sigma_invM) ** 2)

# Regressão

popt, pcov = curve_fit(modelo, invM, y, sigma=erroFinal, absolute_sigma=True)
a_opt, b_opt = popt
perr = np.sqrt(np.diag(pcov))  # erro nos coeficientes
print("erros dos coeficientes a, b:", perr)
print('valores a, b:', a_opt, b_opt)

x_modelo = np.linspace(min(invM), max(invM), 2 * n)
y_modelo = modelo(x_modelo, a_opt, b_opt)

# primeiro gráfico
plt.ylabel(r"$2\Delta S/t^2 \;\, (m/s^2)$")
plt.xlabel(r"$1/m \;\, (kg^{-1})$")
plt.plot(invM, y, 'o', markersize=5, color="black", label="Dados Originais")
plt.plot(x_modelo, y_modelo, color="red", label="Regressão Linear")
plt.errorbar(invM, y, yerr=erroFinal, fmt='o', color="black", capsize=5, label="Erros")
plt.legend(loc="upper left")
plt.show()

# segundo Gráfico, resíduos
res = y - modelo(invM, a_opt, b_opt)

plt.ylabel(r"resíduos $(m/s^2)$")
plt.xlabel(r"$1/m \;\, (kg^{-1})$")
plt.scatter(invM, res, color='b')
plt.errorbar(invM, res, yerr=erroFinal, fmt='o', color="black", capsize=5, label="Erros")
plt.axhline(y=0, color='r', linestyle='--')
plt.show()

# Hipótese nula todos os dados
chi2 = np.sum((y - modelo(invM, a_opt, b_opt)) ** 2 / erroFinal ** 2)
g = len(y) - 2
chi2_red = chi2 / g
print('chi^2 =', chi2_red)

p_chi2 = stats.distributions.chi2.sf(chi2, g)
print("P H_0: ", p_chi2 * 100, " %\n")

# Análise com os gráficos individuais:

invM_1 = []
y_1 = []
y_2 = []
errof_1 = []
errof_2 = []
sigma_y_1 = []
sigma_y_2 = []
sigma_invM_1 = []
sigma_invM_2 = []

for i in range(0, 2 * n, 2):
    invM_1.append(invM[i])
    y_1.append(y[i])
    y_2.append(y[i + 1])
    sigma_y_1.append(sigma_y[i])
    sigma_y_2.append(sigma_y[i + 1])
    sigma_invM_1.append(sigma_invM[i])
    sigma_invM_2.append(sigma_invM[i + 1])

invM_2 = invM_1

# DADOS S1

popti1, pcovi1 = curve_fit(modelo, invM_1, y_1, sigma=sigma_y_1, absolute_sigma=True)
a_i1, b_i1 = popti1

sigma_y_1 = np.array(sigma_y_1)
sigma_invM_1 = np.array(sigma_invM_1)
print('a para o calculo do ErroFinal (só deslocamento 1)', round(a_i1, 2))
errof_1 = np.sqrt((sigma_y_1) ** 2 + a_i1 * a_i1 * (sigma_invM_1) ** 2)
errof_1 = np.array(errof_1)

popt1, pcov1 = curve_fit(modelo, invM_1, y_1, sigma=errof_1, absolute_sigma=True)
a_1, b_1 = popt1
perr1 = np.sqrt(np.diag(pcov1))  # erro nos coeficientes

print("erros a, b (só deslocamento 1):", perr1)
print('a, b (só deslocamento 1)', a_1, b_1)

x_modelo1 = np.linspace(min(np.array(invM_1)), max(np.array(invM_1)), n)
y_modelo1 = modelo(x_modelo1, a_1, b_1)

# gráfico dados s1
plt.ylabel(r"$2\Delta S/t^2 \;\, (m/s^2)$")
plt.xlabel(r"$1/m \;\, (kg^{-1})$")
plt.plot(invM_1, y_1, 'o', markersize=5, color="black", label=r"Dados Originais $\Delta S_1$")
plt.plot(x_modelo1, y_modelo1, color="red", label="Regressão Linear")
plt.errorbar(invM_1, y_1, yerr=errof_1, fmt='o', color="black", capsize=5, label="Erros")
plt.legend(loc="upper left")
plt.show()

# Hipótese nula para s1
chi2_1 = np.sum((y_1 - modelo(np.array(invM_1), a_1, b_1)) ** 2 / errof_1 ** 2)
g1 = len(y_1) - 2
chi2_red_1 = chi2_1 / g1
print('chi^2 reduzido =', chi2_red_1)
p_chi2_1 = stats.distributions.chi2.sf(chi2_1, g1)
print("P H_0 (só deslocamento 1): ", p_chi2_1 * 100, " %\n")

# DADOS S2

popti2, pcovi2 = curve_fit(modelo, invM_2, y_2, sigma=sigma_y_1, absolute_sigma=True)
a_i2, b_i2 = popti2

sigma_y_2 = np.array(sigma_y_2)
sigma_invM_2 = np.array(sigma_invM_2)
print('a para o calculo do ErroFinal (só deslocamento 2)', a_i2)
errof_2 = np.array(np.sqrt(sigma_y_2 ** 2 + a_i2 * a_i2 * sigma_invM_2 ** 2))

popt2, pcov2 = curve_fit(modelo, invM_2, y_2, sigma=sigma_y_2, absolute_sigma=True)
a_2, b_2 = popt2
perr2 = np.sqrt(np.diag(pcov2))  # erro nos coeficientes
print("erros dos coeficientes a, b (só deslocamento 2):", perr2)
print('a, b (só deslocamento 2):', a_2, b_2)

x_modelo2 = np.linspace(min(np.array(invM_2)), max(np.array(invM_2)), n)
y_modelo2 = modelo(x_modelo2, a_2, b_2)

# gráfico dados s2
plt.ylabel(r"$2\Delta S/t^2 \;\, (m/s^2)$")
plt.xlabel(r"$1/m \;\, (kg^{-1})$")
plt.plot(invM_2, y_2, 'o', markersize=5, color="black", label=r"Dados Originais $\Delta S_2$")
plt.plot(x_modelo2, y_modelo2, color="red", label="Regressão Linear")
plt.errorbar(invM_2, y_2, yerr=errof_2, fmt='o', color="black", capsize=5, label="Erros")
plt.legend(loc="upper left")
plt.show()

# Hipótese nula dados S2
chi2_2 = np.sum((y_2 - modelo(np.array(invM_2), a_2, b_2)) ** 2 / errof_2 ** 2)
g2 = len(y_2) - 2
chi2_red_2 = chi2_2 / g2
print('chi^2 =', chi2_red_2)

p_chi2_2 = stats.distributions.chi2.sf(chi2_2, g2)
print("P H_0 (só o deslocamento 2): ", p_chi2_2 * 100, " %\n")

# printar os dados uteis para o relatorio
# print('tempos medios')
# print(', '.join([str(x) for x in tempos]))
#
# print(', '.join([str(x) for x in y]))
#
# print('erros finais')
# print(', '.join([str(x) for x in erroFinal]))
#
# print('y')
# print(', '.join([str(x) for x in y]))
#
# print('y')
# print(', '.join([str(x) for x in y]))
#
# print('inversos das massas')
# print(', '.join([str(x) for x in invM]))
#
# print('erros aceleracoes')
# print(', '.join([str(x) for x in sigma_y]))


# modelo quadratico para a regressao, que é mais verossímil com o experimento,
# uma vez que, sempre há um impulso inicial na soltura do carrinho.

def modeloReal(z, g, v):
    return z * z * g / (2 * (1 / (0.051 * invM) + 1)) + v * z


poptR, pcovR = curve_fit(modeloReal, tempos, deslocamentos, sigma=0.0005)
a_optR, b_optR = poptR
perrR = np.sqrt(np.diag(pcovR))  # erro nos coeficientes
print("erros gReal, vReal:", perrR)
print('gReal, vReal:', a_optR, b_optR)

chi2R = np.sum((deslocamentos - modeloReal(tempos, a_optR, b_optR)) ** 2 / 0.02 ** 2)
gR = len(y) - 2
chi2_redR = chi2R / gR
print('chi^2 =', chi2_redR)

p_chi2R = stats.distributions.chi2.sf(chi2R, gR)
print("P H_0: ", p_chi2R * 100, " %\n")
