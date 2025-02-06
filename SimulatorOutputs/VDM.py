import math
import numpy as np
import whlPhy
import pandas as pd
from matplotlib import pyplot as plt

global tf, tr, lf, lr, L, S, TF, TR, rtbp5
rtbp5 = 0
L = 3
TF = 1.4
TR = 1.4
tf = 0.7
tr = 0.7
lf = 1.4
lr = 1.6
I = 2000


def load_data(row):
    
    v = [row['VX'], row['VY']]
    psi = row['R']
    omega = [row['omega_1'], row['omega_2'], row['omega_3'], row['omega_4']]
    steering = row['STRA']
    lamda = [row['slip_1'], row['slip_2'], row['slip_3'], row['slip_4']]
    tau = [row['TAU_1'], row['TAU_2'], row['TAU_3'], row['TAU_4']]
    # FZ = [row['FFL'],row['FFR'],row['FRL'],row['FRR']]
    # Fext = [row['FZFL'],row['FZFR'],row['FZRL'],row['FZRR']]
    rho = [row['Fz_1'], row['Fz_2'], row['Fz_3'], row['Fz_4']]
    return v, psi, omega, steering, lamda, tau, rho

def force(FFL, FFR, FRL, FRR, SA): 

    dl = SA[0]
    dr = SA[1]

    fx = FFL[0]*math.cos(dl) + FFL[1]*math.sin(dl) + FFR[0]*math.cos(dr) + FFR[1]*math.sin(dr) + FRL[0] + FRR[0]
    fy = FFL[0]*math.sin(dl) - FFL[1]*math.cos(dl) + FFR[0]*math.sin(dr) - FFR[1]*math.cos(dr) - FRL[1] - FRR[1]

    return [fx, fy]

def moment(FFL, FFR, FRL, FRR, SA):

    dl = SA[0]
    dr = SA[1]

    xm = tf*(FFL[0]*math.cos(dl) - FFL[1]*math.sin(dl)) - tf*(FFR[0]*math.cos(dr) - FFR[1]*math.sin(dr)) + tr*FRL[0] - tr*FRR[0]
    ym = lf*(FFL[0]*math.sin(dl) + FFL[1]*math.cos(dl)) + lf*(FFR[0]*math.sin(dr) + FFR[1]*math.cos(dr)) - lr*FRL[1] - lr*FRR[1]

    m = xm + ym

    return m


def create_csv(ax,ay,ag,ti):
    data = {
    'TIME': ti,    
    'AX_NEW': ax,
    'AY_NEW': ay,
    'AngACC_NEW': ag
    }

    # Create a DataFrame
    df = pd.DataFrame(data)

    # Specify the CSV file name
    csv_file = 'prediction.csv'

    # Write DataFrame to CSV
    df.to_csv(csv_file, index=False)

    print(f'CSV file "{csv_file}" has been created successfully.')

def create_csv2(ti,LSFL,LSFR,LSRL,LSRR,SAFL,SAFR,SARL,SARR):
    data = {
    'TIME': ti,    
    'lamdaFL': LSFL,
    'lamdaFR': LSFR,
    'lamdaRL': LSRL,
    'lamdaRR': LSRR,
    'fxFL': SAFL,
    'tauFR': SAFR,
    'tauRL': SARL,
    'tauRR': SARR,
    }

    # Create a DataFrame
    df = pd.DataFrame(data)

    # Specify the CSV file name
    csv_file = 'training.csv'

    # Write DataFrame to CSV
    df.to_csv(csv_file, index=False)

    print(f'CSV file "{csv_file}" has been created successfully.')


def plot_lamda(lamda,lamda2):
    fig, axs = plt.subplots(2,1)
    axs[0].plot(lamda, label = 'actual')
    axs[0].plot(lamda2, label = 'prediction')
    # axs[1].plot(om, label = 'omegaX')
    axs[0].legend()
    # axs[1].legend()
    # axs[0].plot(axx, label='acceleration')
    plt.show()

def calculate_fext(ax,ay,v,r,FFL,FFR,FRL,FRR):

    b = 1.6
    a = 1.4
    m = 1500
    g = 9.81
    h = 0.35
    wf = 1.4
    d = 0
    vx = v[0]
    vy = v[1]

    fzf = ((b*m*g)/(a+b)) -  (((ax)*m*h*9.81)/(a + b))
    fzr = ((a*m*g)/(a+b)) +  (((ax)*m*h*9.81)/(a + b))

    fzfl = fzf/2 + (m*h*(ay))*4/wf 
    fzfr = fzf/2 + (-m*h*(ay))*4/wf
    fzrl = fzr/2 + (m*h*(ay))*4/wf  + 90*9.81
    fzrr = fzr/2 + (-m*h*(ay))*4/wf + 90*9.81

    # fext = [fzfl,fzfr,fzrl,fzrr]


    # rtbfi5 = v[1]*r
    # rtbvcp = m*g
    # rtbfi8 = 0
    # rtbp4i1 = 0
    # rtbp7 = 0 #1?
    # rtbfi11 = v[0]*r
    # bb1tmp = ay
    # bb1tmp_0 = bb1tmp*m  #(m*ay)
    # bB1 = rtbp7 - (0 - bb1tmp_0)*h  #(m*ay*h)
    # rtbp5ttmp = ax 
    # rtbp5temp0 = rtbp5ttmp*m #(m*ax)
    # rtbp5 = (((0 - rtbp5temp0) - rtbfi8)*h - 0*h) - rtbp4i1 #(-m*ax*h)
    # bb3tmp = -rtbvcp #(-m*g)
    # rtbfi2 = a + b

    # b1 = FFL[0]
    # b2 = FFR[1]
    # b3 = FRL[0]

    # Fztmp = bb3tmp*b #(-m*g*b)
    # fztmp0 = 2*bB1*a + 2*bB1*b #(m*ay*h)*(a+b)*2
    # fztmp1 = rtbp5*w # (-m*ax*h)*w
    # fztmp2 = rtbp5*w # (-m*ax*h)*w
    # fztmp3 = 2*bb3tmp*a*d #-2*m*g*a*d or 0
    # fztmp4 = 2*bb3tmp*b*d
    # fztmp5 = Fztmp*w #(-m*g*b*w)
    # Fztmp *= w #(-m*g*b*w)
    # fztmp6 = rtbfi2*2*(w+w) #4*W*(a+b))
    # fztmp7 = (((fztmp0 + fztmp1) + fztmp2) + fztmp3) + fztmp4  # (m*ay*h)*(a+b)*2 + (-m*ax*h)*w + (-m*ax*h)*w + 0 + 0 
    # fzfl = -(-((fztmp7 - fztmp5) - Fztmp)/fztmp6)*2 # -(-(2((m*ay*h)*(a+b) - m*ax*h*w) - (-m*g*b*w) - (-m*g*b*w))/((-m*g*b*w)))*2
    # fztmp0 = (((fztmp0 - fztmp1) - fztmp2) + fztmp3) + fztmp4
    # fzfr = -(((fztmp0 + fztmp5) + Fztmp)/fztmp6)*2
    # Fztmp = bb3tmp*a
    # fztmp1 = Fztmp*w
    # Fztmp *= w
    # fzrl = -(-((fztmp0 - fztmp1) - Fztmp)/fztmp6)*2
    # fzrr = -(((fztmp7 + fztmp1) + Fztmp)/fztmp6)*2

    fext = [fzfl,fzfr,fzrl,fzrr]
    
    for i in range(len(fext)):
        if fext[i] < 0:
            fext[i] = 0
        
        i += 1

    return fext
      
input = pd.read_csv('input_common_sample.csv')
ro = 0.30454
rl = ro/1.00866439868088


# wheelFL = whlPhy.Wheel(0.322361262,'FL')
# wheelFL = whlPhy.Wheel(0.32215,'FL')
# wheelFR = whlPhy.Wheel(0.32248113,'FR')
# wheelRL = whlPhy.Wheel(0.31966804,'RL')
# wheelRR = whlPhy.Wheel(0.319347774,'RR')

wheelFL = whlPhy.Wheel(ro,'FL')
wheelFR = whlPhy.Wheel(ro,'FR')
wheelRL = whlPhy.Wheel(ro,'RL')
wheelRR = whlPhy.Wheel(ro,'RR')

ax = []
ay = []
ti = []
longslipFL = []
longslipFR = []
longslipRL = []
longslipRR = []
fxFL = []
tauFR = []
tauRL = []
tauRR = []
alpha = []
FFL = [0,0]
FFR = [0,0]
FRL = [0,0]
FRR = [0,0]
AX = 0
AY = 0

for index,row in input.iterrows():
    v, psi, omega, steering, lamda, tau, r = load_data(row)

    fe = calculate_fext(AX/9.8,AY/9.8,v,psi,FFL,FFR,FRL,FRR)

    wheelFL.update(v, psi, lamda[0], steering, omega[0],tau[0],row['TIME'],r[0])
    wheelFR.update(v, psi, lamda[1], steering, omega[1],tau[1],row['TIME'],r[1])
    wheelRL.update(v, psi, lamda[2], steering, omega[2],tau[2],row['TIME'],r[2])
    wheelRR.update(v, psi, lamda[3], steering, omega[3],tau[3],row['TIME'],r[3])

    longslipFL.append(wheelFL._Wheel__lamda)
    longslipFR.append(wheelFR._Wheel__lamda)
    longslipRL.append(wheelRL._Wheel__lamda)
    longslipRR.append(wheelRL._Wheel__lamda)
    tauFR.append(wheelFR._Wheel__tau)
    tauRL.append(wheelRL._Wheel__tau)
    tauRR.append(wheelRR._Wheel__tau)

    FFL = wheelFL.Pacejka()
    FFR = wheelFR.Pacejka()
    FRL = wheelRL.Pacejka()
    FRR = wheelRR.Pacejka()
# 
    # FFL = [row['FXFL'],row['FYFL']]
    # FFR = [row['FXFR'],row['FYFR']]
    # FRL = [row['FXRL'],row['FYRL']]
    # FRR = [row['FXRR'],row['FYRR']]

    # fxFL.append((fe[0]*1.93553586)-1.12672731E4)
    fxFL.append(wheelFL._Wheel__Fz)


    steeringA = wheelFL.steering_angle(steering,both=True)
    
    forces = force(FFL, FFR, FRL, FRR, steeringA)
    # forces = [(row['FXFL']+row['FXFR']+row['FXRL']+row['FXRR']),(row['FYFL']+row['FYFR']+row['FYRL']+row['FYRR'])]
    Zmoment = moment(FFL, FFR, FRL, FRR, steeringA)

    AX = forces[0]/1500
    AY = forces[1]/1500
    Za = Zmoment/I
    
    ax.append(((AX)/9.8)-0.01)
    ay.append(((AY)/9.8)+0.025)
    alpha.append(-Za)
    ti.append(row['TIME'])


create_csv(ax,ay,alpha,ti)
create_csv2(ti,longslipFL,longslipFR,longslipRL,longslipRR,fxFL,tauFR,tauRL,tauRR)
# plot_lamda(input['FFL'] - input['FZFL'],0)

    

