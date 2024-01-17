import numpy as np

def calcular(a, b, c, d, f, engranaje, phi, theta_2, omega_2, alfa_2, p, delta_3):

    #CÁLCULO DE theta_5
    theta_5 = engranaje*theta_2 + phi

    #VARIABLES
    A = 2*c*(d*np.cos(np.radians(theta_5))-a*np.cos(np.radians(theta_2))+f)
    B = 2*c*(d*np.sin(np.radians(theta_5))-a*np.sin(np.radians(theta_2)))
    C = a**2-b**2+c**2+d**2+f**2-2*a*f*np.cos(np.radians(theta_2))-2*d*(a*np.cos(np.radians(theta_2))-f)*np.cos(np.radians(theta_5))-2*a*d*np.sin(np.radians(theta_2))*np.sin(np.radians(theta_5))
    D = C-A
    E = 2*B
    F = A+C

    #CÁLCULO DE theta_4
    theta_4_1 = 2*np.degrees(np.arctan((-E+np.sqrt(E**2-4*D*F))/(2*D)))
    theta_4_2 = 2*np.degrees(np.arctan((-E-np.sqrt(E**2-4*D*F))/(2*D)))

    #CÁLCULO DE theta_3
    theta_3_1 = np.arcsin((-a*np.sin(np.radians(theta_2))+c*np.sin(np.radians(theta_4_1))+d*np.sin(np.radians(theta_5)))/b)
    theta_3_2 = np.arcsin((-a*np.sin(np.radians(theta_2))+c*np.sin(np.radians(theta_4_2))+d*np.sin(np.radians(theta_5)))/b)

    #CÁLCULO DE omega_5
    omega_5 = engranaje*omega_2

    #CÁLCULO DE omega_3 Y omega_4
    theta_3 = theta_3_1
    theta_4 = theta_4_1

    omega_3 = -(2*np.sin(np.radians(theta_4))*((a*omega_2*np.sin(np.radians(theta_2-theta_4))+d*omega_5*np.sin(np.radians(theta_4-theta_5)))))/(b*(np.cos(np.radians(theta_3-2*theta_4))-np.cos(np.radians(theta_3))))
    omega_4 = (a*omega_2*np.sin(np.radians(theta_2))+b*omega_3*np.sin(np.radians(theta_3))-d*omega_5*np.sin(np.radians(theta_5)))/(c*np.sin(np.radians(theta_4)))

    #CÁLCULO DE VELOCIDADES
    V_Ax = -a*omega_2*np.sin(np.radians(theta_2))
    V_Ay = a*omega_2*np.cos(np.radians(theta_2))
    V_A = np.sqrt((V_Ax)**2+(V_Ay)**2)

    V_BAx = -b*omega_3*np.sin(np.radians(theta_3))
    V_BAy = b*omega_3*np.cos(np.radians(theta_3))
    V_BA = np.sqrt((V_BAx)**2+(V_BAy)**2)

    V_Bx = V_Ax + V_BAx
    V_By = V_Ay + V_BAy
    V_B = np.sqrt((V_Bx)**2+(V_By)**2)

    V_Cx = -d*omega_5*np.sin(np.radians(theta_5))
    V_Cy = d*omega_5*np.cos(np.radians(theta_5))
    V_C = np.sqrt((V_Cx)**2+(V_Cy)**2)

    V_PAx = -p*omega_3*np.sin(np.radians(theta_3+delta_3))
    V_PAy = p*omega_3*np.cos(np.radians(theta_3+delta_3))
    V_PA = np.sqrt((V_PAx)**2+(V_PAy)**2)

    V_Px = V_Ax + V_PAx
    V_Py = V_Ay + V_PAy
    V_P = np.sqrt((V_Px)**2+(V_Py)**2)

    #CÁLCULO DE alfa_5
    alfa_5 = engranaje*alfa_2

    #CÁLCULO DE alfa_3 Y alfa_4
    alfa_3 = (-a*alfa_2*np.sin(np.radians(theta_2-theta_4))-a*(omega_2)**2*np.cos(np.radians(theta_2-theta_4))-b*(omega_3)**2*np.cos(np.radians(theta_3-theta_4))+d*(omega_5)**2*np.cos(np.radians(theta_5-theta_4))+d*alfa_5*np.sin(np.radians(theta_5-theta_4))+c*(omega_4)**2)/(b*np.sin(np.radians(theta_3-theta_4)))
    alfa_4 = (a*alfa_2*np.sin(np.radians(theta_2-theta_3))+a*(omega_2)**2*np.cos(np.radians(theta_2-theta_3))-c*(omega_4)**2*np.cos(np.radians(theta_3-theta_4))-d*(omega_5)**2*np.cos(np.radians(theta_3-theta_5))+d*alfa_5*np.sin(np.radians(theta_3-theta_5))+b*(omega_3)**2)/(c*np.sin(np.radians(theta_4-theta_3)))

    #CÁLCULO DE ACELERACIONES
    A_Ax = -a*alfa_2*np.sin(np.radians(theta_2))-a*(omega_2)**2*np.cos(np.radians(theta_2))
    A_Ay = a*alfa_2*np.cos(np.radians(theta_2))-a*(omega_2)**2*np.sin(np.radians(theta_2))
    A_A = np.sqrt((A_Ax)**2+(A_Ay)**2)

    A_BAx = -b*alfa_3*np.sin(np.radians(theta_3))-b*(omega_3)**2*np.cos(np.radians(theta_3))
    A_BAy = b*alfa_3*np.cos(np.radians(theta_3))-b*(omega_3)**2*np.sin(np.radians(theta_3))
    A_BA = np.sqrt((A_BAx)**2+(A_BAy)**2)

    A_Bx = A_Ax + A_BAx
    A_By = A_Ay + A_BAy
    A_B = np.sqrt((A_Bx)**2+(A_By)**2)

    A_Cx = -c*alfa_5*np.sin(np.radians(theta_5))-c*(omega_5)**2*np.cos(np.radians(theta_5))
    A_Cy = c*alfa_5*np.cos(np.radians(theta_5))-c*(omega_5)**2*np.sin(np.radians(theta_5))
    A_C = np.sqrt((A_Cx)**2+(A_Cy)**2)

    A_PAx = -p*alfa_3*np.sin(np.radians(theta_3+delta_3))-p*(omega_3)**2*np.cos(np.radians(theta_3+delta_3))
    A_PAy = p*alfa_3*np.cos(np.radians(theta_3+delta_3))-p*(omega_3)**2*np.sin(np.radians(theta_3+delta_3))
    A_PA = np.sqrt((A_PAx)**2+(A_PAy)**2)

    A_Px = A_Ax + A_PAx
    A_Py = A_Ay + A_PAy
    A_P = np.sqrt((A_Px)**2+(A_Py)**2)

    resultados = [alfa_3, alfa_4, A_A, A_B, A_C, A_P]
    print(resultados)

fila_a = [1, 7, 9, 4, 6, 2.0, 30, 60, 10, 0, 6, 30]
fila_b = [5, 7, 8, 4, 6, -2.5, 60, 30, -12, 5, 9, 25]
fila_c = [5, 7, 8, 4, 3, -0.5, 0, 45, -15, -10, 10, 80]
fila_d = [5, 7, 8, 4, 4, -1.0, 120, 75, 24, -4, 5, 45]
fila_e = [9, 11, 8, 8, 5, 3.2, -50, -39, -50, 10, 9, 300]
fila_f = [2, 7, 5, 3, 10, 1.5, 30, 120, -45, 50, 10, 120]
fila_g = [7, 9, 11, 4, 15, 2.5, -90, 75, 100, 18, 4, 300]
fila_h = [8, 7, 9, 4, 12, -2.5, 60, 55, -65, 25, 6, 20]
fila_i = [7, 8, 9, 4, 9, -4.0, 120, 100, 25, -25, 9, 80]


calcular(1, 7, 9, 4, 6, 2.0, 30, 60, 10, 0, 6, 30)
calcular(7, 8, 9, 4, 9, -4.0, 120, 100, 25, -25, 9, 80)
