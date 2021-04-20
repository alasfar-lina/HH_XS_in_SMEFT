# -*- coding: utf-8 -*-
import rundec  # required for alpha_s running
import numpy as np
crd = rundec.CRunDec()
## alpha_s running stuff
alps_mz = 0.1183
Mz = 91.1876
nf = 5

def funcalps(muR):
  return  crd.AlphasExact(alps_mz, Mz,muR, nf, 3)

# All masses in GeV.

v=246.0
MT= 173.2
MH=125.1
#alps= 0.1198 # alpha_s at MZ, it is better to use the running alphas
alps =  funcalps(2.*MH) # renormalisation scale. 
L=1e6 # new physics scale


def XSNLO(CH,CgH,ckin,CtH):
  cH=CH/L
  cuH=CtH/L
  cHkin=ckin/L
  cHG=CgH/L/16/np.pi**2
  return 1.000-3.9866459999999995*cHkin*v**2+8.949195999999999*cHkin**2*v**4 +((3.1233253007332573 -14.981140447025783*cHkin*v**2)*cuH*v**3)/MT+(cH*(1.7820680000000002-4.268759999999999*cHkin*v**2)*v**4)/MH**2 +(6.637391499999998*cuH**2*v**6)/MT**2 +(2.545778159529616*cH*cuH*v**7)/(MH**2*MT) + (1.368992*cH**2*v**8)/MH**4+ (3106.4364774009646*cHG**2*v**4)/alps**2 - (74.20004538081712*cH*cHG*v**6)/(MH**2*alps)+(cHG*v**2*((-52.47197151255347 + 236.21675299329388*cHkin)*MT -188.14684359273673*cuH*v**3))/(MT*alps)



def XSLO(CH,CgH,ckin,CtH):
  cH=CH/L
  cuH=CtH/L
  cHkin=ckin/L
  cHG=CgH/L/16/np.pi**2
  return 1.-3.2954600000000003*cHkin*v**2 +6.722859999999997*cHkin**2*v**4+((2.6399053858140427- 11.480309927700123*cHkin*v**2)*cuH*v**3)/MT+(cH*(1.6048999999999998-3.015439999999998*cHkin*v**2)*v**4)/MH**2+(5.18126*cuH**2*v**6)/MT**2+(1.753815736173558*cH*cuH*v**7)/(MH**2*MT)+(1.11256*cH**2*v**8)/MH**4+(3050.5592304475026*cHG**2*v**4)/alps**2-(60.40684513611929*cH*cHG*v**6)/(MH**2*alps)+(cHG*v**2*((-53.14725283220398+186.0102075680209*cHkin)*MT-154.7657478815316*cuH*v**3))/(MT*alps)


def Kfac(CH,CgH,ckin,CtH):
  num =XSNLO(CH,CgH,ckin,CtH)
  den=XSLO(CH,CgH,ckin,CtH)
  return (num/den)*1.6609588715365242




if __name__ == '__main__':
    print("The k-factor for CH= 10 is %f"%(Kfac(10.0,0.0,0.0,0.0)))
