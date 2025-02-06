import numpy as np
import math
from scipy.integrate import cumtrapz

class Wheel:
    
    global tf, lf, B, C, D, E, L, S, TF, TR, tr, tf, lr,Sv
    S = 6
    L = 3
    TF = 1.4
    TR = 1.4
    tf = 0.7
    tr = 0.7
    lf = 1.4
    lr = 1.6
    Sv = 230

    def __init__(self,radius, location):
        self.__radius = radius
        self.__location = location
        self.__VR = None
        self.__lamda = None
        self.__lamda2 = None
        self.__tau = None
        self.__velprev = 0
        self.__vel = np.array([0])
        self.__val = True
        self.__Fx = 0
        self.__Fy = 0
        if location == 'FL':
            self.__mass = 350*9.8
        elif location == 'FR':
            self.__mass = 350.0*9.8
        elif location == 'RL':
            self.__mass = 400.0*9.8
        else:
            self.__mass = 400.0*9.8

        self.__Fz = self.__mass
 
    def steering_angle(self, steering, both):
        steering = steering + 1E-10
        if both == True:
            delL = math.atan(L/((L/(math.tan(steering/S))) - (TF/2)))
            delR = math.atan(L/((L/(math.tan(steering/S))) + (TF/2)))
            str = [delL,delR]
            return str
        else:
            if self.__location == 'FL':
                str = math.atan(L/((L/(math.tan(steering/S))) - (TF/2)))
                return str
            elif self.__location == 'FR':
                str = math.atan(L/((L/(math.tan(steering/S))) + (TF/2)))
                return str
            else:
                print("This is a rear tyre, aborting calculation!")
                return None

    def rolling_radius(self):

        qre0 = 1.00866439868088
        qv1 = 0.000760413786224011
        fnom = 4600
        Dreff = 0.260468730454265
        Breff = 8.25094594147963
        Freff = 0.0735298544471851
        pfz1 = 0.69958166705601
        pi = 2.4e5
        pnom = 262000
        dpi = (pi - pnom)/pnom
        qfcy2 = 0.10843499565426
        ro = 0.33136
        longvl = 16.7
        qfz1 = 0
        qfz2 = 15.6870832810226

        # cz = (qfcy2*self.__Fy/fnom + self.__Fz/self.__rho)
        # cz = (fnom/ro)*(qfz1 + 2*qfz2*(self.__rho/ro))
        # cz *= (pfz1*dpi + 1)
        # self.__cz = cz
        cz = 172408.4927

        dc_data = self.__omega*ro/longvl
        lcz = self.__Fz/fnom

        re = ((dc_data**2)*qv1 + qre0)*ro -(math.atan(Breff*lcz)*Dreff + Freff*lcz)*(fnom/cz)

        return re

            
    def rotate(self,ang):
        x = self.__VR[0]*math.cos(ang) + self.__VR[1]*math.sin(ang)
        y = self.__VR[1]*math.cos(ang) - self.__VR[0]*math.sin(ang)

        return [x,y]
    
    def find_params_x(self):

        Fzo = 4600
        Fn = 4600
        Fz = self.__Fz
        dfz = (Fz - Fzo)/Fzo
        pi = 2.4e5
        pn = 262000
        dpi = (pi - pn)/pn
        # dfz = 0
        self.__dfz = dfz
        vo = math.sqrt(9.8*0.33136)
        ro = 0.33136
        rl = ro/1.00866439868088
        qv1 = 0.000760413786224011
        qv2 = 0.0463384792019201
        qfz1 =  0
        qfz2 = 15.6870832810226
        pvx1 = 1.59236164471432e-05
        pvx2 = 0.000104321112127671
        phx1 = 0.000212599305364818
        phx2 = 0.00115950515263055
        pkx1 = 21.3872544684023
        pkx2 = 14.0006541873175
        pkx3 = -0.405326109653452
        pex1 = 0.110312777158455
        pex2 = 0.313468516038146
        pex3 = 0
        pex4 = 0.0016060905760045
        pdx1 = 1.0239116238178
        pdx2 = -0.0842405110022724
        pcx1 = 1.58523057950359
        qfcx1 = 0.138643970247602
        qfcy1 = 0.10843499565426
        ppx1 = -0.349461321276586
        ppx2 = 0.387840040616429
        ppx3 = -0.0969947336569324
        ppx4 = 0.0632271859795801
        ev = 0.1

        dr = (qv1*ro)*(((self.__omega*ro)/vo)**2)
        rc = (ro - rl + dr)
        rhoz = max(rc,0)

        Svx = Fz*(pvx1 + pvx2*dfz)

        Shx = (phx1 + phx2*dfz)

        kx = self.__lamda + Shx

        Kxk = Fz*(pkx1 + pkx2*dfz)*np.exp(pkx3*dfz)*(1 + ppx1*dpi + ppx2*(dpi**2))

        E = (pex1 + pex2*dfz +pex3*dfz**2)*(1-pex4*np.sign(kx))

        ux = (pdx1 + pdx2*dfz)*(1 + ppx3*dpi + ppx4*(dpi**2))

        D = ux*Fz
        C = pcx1
        C = max(C,0)
        D = max(D,0)
        E = min(E,1)
        B = Kxk/(C*D + ev)

        return B,C,D,E,Svx,kx
    
    def find_params_y(self):

        pi = 2.4e5
        pn = 262000
        dpi = (pi - pn)/pn
        vo = math.sqrt(9.8*0.33136)
        ro = 0.33136
        rl = ro/1.00866439868088
        Fzo = 4600
        pky1 = -15.5714726518315
        pky2 = 1.73126522291751
        pky4 = 1.98176755955416
        pky6 = -0.884005199550313
        pky7 = -0.237259727611847
        pvy1 = -0.00675427560163264
        pvy2 = 0.036379218103843
        phy1 = -0.00183370557235628
        phy2 = 0.00346401302716671
        ek = 0.01
        pdy1 = 0.878267729082195
        pdy2 = -0.0644597923147385
        pcy1 = 1.34299950037631
        ey = 0.1
        pey1 = -0.809776534470972
        pey2 = -0.600180598867198
        pey3 = 0.0991732552633487
        ppy1 = -0.62059646553185
        ppy2 = -0.0647824842338686
        ppy3 = -0.164648843290686
        ppy4 = 0.283193909060098
        Fz = self.__Fz 
        dfz = (Fz - Fzo)/Fzo
    

        # dfz = 0

        Svyg = 0
        Kya = pky1*Fzo*(1 + ppy1*dpi)*math.sin(pky4*math.atan(Fz/(pky2*Fzo*(1 + ppy2*dpi))))
        Kygo = Fz*(pky6 + pky7*dfz)
        Svy = Fz*(pvy1 + pvy2*dfz) 
        Shy = (phy1 + phy2*dfz) + (-Svyg)/(Kya + ek)
        uy = (pdy1 + pdy2*dfz)*(1 + ppy3*dpi + ppy4*(dpi**2))

        Dy = uy*Fz
        Cy = pcy1
        Cy = max(Cy,0)
        ay = self.__tau + Shy
        By = Kya/(Cy*Dy + ey)
        Ey = (pey1 + pey2*dfz)*(1 - pey3*np.sign(ay))
        Ey = min(Ey,1)

        return By,Cy,Dy,Ey,Svy,ay


    
    def combined(self,fx,fy):

        pi = 2.4e5
        pn = 262000
        rhx1 = 0.000233272373115809
        rex1 = -0.45202516851367
        rex2 = -0.47304886171176
        rcx1 = 1.02796280922059
        rbx1 = 12.7633329850276
        rbx2 = 9.5787123658471
        rhy1 = 0.0118479204960415
        rcy1 = 1.06691874059059
        rey1 = 0.308244924574305
        rby1 = 10.7588809173639
        rby2 = 7.70420065295029
        rby3 = 3.40283755070691e-06
        pey1 = -0.809776534470972
        rvy1 = 0.0560024425050078
        rvy4 = 98.4047026517149
        rvy5 = 2.02759273900054
        rvy6 = 15.7623064923347
        pdy1 = 0.878267729082195
        pdy2 = -0.0644597923147385
        rhy2 = 7.63650177161071e-06
        rey2 = 7.48775032604508e-06
        rvy2 = 7.48487127056197e-06
        ppy1 = -0.62059646553185
        ppy2 = -0.0647824842338686
        ppy3 = -0.164648843290686
        ppy4 = 0.283193909060098

        dpi = (pi - pn)/pn


        Shxa = rhx1

        Exa = rex1 + rex2*self.__dfz

        cxa = rcx1

        Bxa = (rbx1)*math.cos(math.atan(rbx2*self.__lamda))
        Bxa = max(Bxa,0)
        Exa = min(Exa,1)


        As = self.__tau + Shxa
        Gxao = math.cos(cxa*math.atan(Bxa*Shxa-Exa*(Bxa*Shxa-math.atan(Bxa*Shxa))))
        Gxao = max(Gxao,0)
        BA = Bxa*As
        Gxa = (math.cos(cxa*math.atan(BA-Exa*(BA-math.atan(BA)))))/Gxao

        Fx = fx*Gxa

        #########################

        Shyk = rhy1 + rhy2*self.__dfz
        Cyk = rcy1
        Eyk = rey1 + rey2*self.__dfz
        Byk = rby1*(math.cos(math.atan(rby2*(self.__tau - rby3))))
        ks = self.__lamda + Shyk
        BS = Byk*Shyk
        BK = Byk*ks
        uy = (pdy1 + pdy2*self.__dfz)*(1 + ppy3*dpi + ppy4*(dpi**2))
        Dvyk = uy*self.__Fz*(rvy1 + rvy2*self.__dfz)*(math.cos(math.atan(rvy4*(self.__tau))))
        Svyk = Dvyk*math.sin(rvy5*math.atan(rvy6*self.__lamda))

        Byk = max(Byk,0)
        Eyk = min(Eyk,1)
        
        Gyko = math.cos(Cyk*math.atan(BS - Eyk*(BS - math.atan(BS))))
        Gyko = max(Gyko,0)
        Gyk = (math.cos(Cyk*math.atan(BK - Eyk*(BK - math.atan(BK)))))/Gyko

        Fy = Gyk*fy + Svyk


        return Fx,Fy
    

    def rho(self,fext):

        # tot = self.__Fz - 10.35*9.81 - fext - 494.65*self.__velprev

        self.__tot = fext
        

        self.__A = np.append(self.__A,[self.__tot/10.35])
        self.__T = np.append(self.__T,np.array([self.__time]))

        vel = cumtrapz(self.__A,self.__T,0)

        self.__vel = np.append(np.array([0]),vel)
    
        self.__velprev = vel[-1]

        # if self.__velprev < 0:
            # self.__velprev = 0

        rho = cumtrapz(self.__vel,self.__T,0)

        if -rho[-1] < 0:
            rho[-1] = 0

        return -rho[-1]


    def calculate_fz(self):

        qfcx = 0.138643970247602
        qfcy = 0.10843499565426
        fzo = 4600
        qv2 = 0.0463384792019201
        qfz1 = 0
        qfz2 = 15.6870832810226
        pfz1 = 0.69958166705601
        pi = 2.4e5
        pn = 262000
        dpi = (pi - pn)/pn
        # Cfz = 207285.061134007
        Cfz = 183168.446770145
        # Cfz = 172408.4927


        a = ((qfcx*self.__Fx)/fzo)**2
        b = ((qfcy*self.__Fy)/fzo)**2
        ro = 0.33136
        vo = math.sqrt(ro*9.8)
        Fz = max(Cfz*self.__rho,100)

        # Fz = (1 + ((qv2*np.abs(self.__omega)*ro)/vo) - a -b)*(((qfz1*self.__rho)/ro) + ((qfz2*(self.__rho**2))/ro**2))*(1 + pfz1*dpi)*fzo

        return Fz 





            
    def update(self, V, PSI, lamda, steering, omega, tau,time,rho): 
        
        self.__omega = omega

        # self.__Fz = fz 
        self.__time = time
        # self.__rho = self.rho(rho)
        self.__rho = rho
        self.__Fz = self.calculate_fz()
        # self.__Fz = fz

        # if self.__val == True:
            # self.__rho = 0
            # self.__val = False
        # else:
            # self.__rho = self.rho(fext)

        self.__radius2 = self.rolling_radius()

        if self.__location == 'FL':
        
            VX = V[0] + tf*PSI
            VY = V[1] + lf*PSI
        elif self.__location == 'FR':
            VX = V[0] - tf*PSI
            VY = V[1] + lf*PSI  
        elif self.__location == 'RL':
            VX = V[0] + tr*PSI
            VY = V[1] - lr*PSI
        else:
            VX = V[0] - tr*PSI
            VY = V[1] - lr*PSI 
        
        self.__VR = [VX,VY]
    
        if self.__location == 'FL' or self.__location == 'FR':
            ang = self.steering_angle(steering, both=False)
            VT = self.rotate(ang)
            self.__VT = self.rotate(ang)
        else:
            VT = self.__VR
        
# 
        if omega == 0:
            self.__lamda = 0
            self.__tau = 0
        elif math.sqrt(VT[0]**2 + VT[1]**2) <= 0.05:
            self.__lamda = 0
            self.__tau = 0
        else:

        # if math.sqrt(np.abs(VT[0]**2) + np.abs(VT[1]**2)) <= 0.06:
            # self.__lamda = 0
            # self.__tau = 0
        # else:    
            if VT[0] > self.__radius2*omega:
                self.__lamda = ((self.__radius2*omega)-VT[0])/np.absolute(VT[0])
                # self.__lamda = lamda
            else:
                self.__lamda = ((self.__radius2*omega)-VT[0])/np.absolute(VT[0])
                # self.__lamda = lamda

            self.__lamda = min(self.__lamda,1.5)
            self.__lamda = max(self.__lamda,-1.5)
 
            self.__lamda = max(min(self.__lamda,1),-1)
            if self.__lamda > 1 or self.__lamda < -1:
                self.__lamda = lamda
            # self.__lamda = lamda 
#  
        # self.__lamda = lamda
            if self.__location == 'FL' or self.__location == 'FR':
                self.__tau =  ang - math.atan2(self.__VR[1],(self.__VR[0]))
                # self.__tau = tau
            else:
                # self.__tau = tau
                self.__tau = math.atan2(-self.__VR[1],(self.__VR[0]))
            self.__tau = min(self.__tau,1.5)
            self.__tau = max(self.__tau,-1.5)   


        # self.__tau = tau

    def Pacejka(self):

        # D = 0.6
        # B = 20
        # C = 1.9
        # E = 0.97
        # B = 1.0543224400440097
        # C = 8.343198725304742
        # D = 8.868561577644178
        # E = 1.17695912327173
        

        # D = 6.98100915539891s6
        # B = 0.00683218837162s7478
        # C = 93.9688263090983s8
        # E = 11.2344602378345s

    
        # B = 0.011847209814202008
        # C = 30.018269608864624
        # D = 10.981012134995351
        # E = 1.98584724643717

        B = 0.010847209814202008
        C = 30.18269608864624
        D = 10.981012134995351
        E = 1.98584724643717
#  
        # By = 0.010847209814202008
        # Cy = 30.18269608864624
        # Dy = 10.981012134995351
        # Ey = 1.98584724643717

        Bx,Cx,Dx,Ex,Svx,kx = self.find_params_x()
        By,Cy,Dy,Ey,Svy,ay = self.find_params_y()
        # Dx = 1.0537335545532367*self.__Fz
        # Ex = self.__Fz*0.000018885*3.1 - 0.16
        # Bx = self.__Fz**0.00324841144472165*0.35 + 8.398656442
        self.__Bx = Bx
        # self.__Dx = Dx
        # self.__Ex = Ex
        # self.__Cx = Cx
        # self.__svx = Svx
        # Dy = self.__Fz*0.9093912496970589
        # Ey = self.__Fz*-0.00019105152572105878
        # By = self.__Fz*0.0031331122512320074*0.3 - 24 + 8
        # self.__By = By
        # self.__Dy = Dy
        # self.__Ey = Ey
        # self.__Cy = Cy

    
    
     
        BL = Bx*kx
        BT = By*ay
        # BL = Bx*self.__lamda
        # BT = By*(self.__tau)
        
# 
        FTx = Dx*math.sin(Cx*math.atan(BL-Ex*(BL-math.atan(BL)))) + Svx
        FTy = Dy*math.sin(Cy*math.atan(BT-Ey*(BT-math.atan(BT)))) + Svy
# 
        self.__Fx,self.__Fy = self.combined(FTx,FTy)

        # FTx = self.__Fz*D*math.sin(C*math.atan(BL-E*(BL-math.atan(BL))))
        # FTy = self.__Fz*D*math.sin(C*math.atan(BT-E*(BT-math.atan(BT))))

        return [self.__Fx,self.__Fy]


        
    

    
    




