import yaml
filename ="./data.yaml"
global stream
stream = open(filename, 'r')

# calculate the branching ratio for bbx aa as a function of kappa_q
# f= 'ku','kd'..etc
def BR(kq,f):
    data= yaml.safe_load(stream)
    Rga= data[f]['Width']['gagaA0']+ data[f]['Width']['gagaA1']*kq+data[f]['Width']['gagaA2']*kq**2
    Rtot= data[f]['Width']['totA0']+ data[f]['Width']['totA1']*kq+data[f]['Width']['totA2']*kq**2
    HiggsFullwidth = data['Higgs']['width']*Rtot
    return 2*(data['Higgs']['BRbbSM']*Rga*data['Higgs']['BRgagaSM'])/Rtot


# calculate the NLO cross-section for qq -> hh as a function of kappa_q
#W ='13TeV', '14TeV', '27TeV','100TeV'
def XS(kq,f,W):
    data= yaml.safe_load(stream)
    R=data[f][W]['XSA0']+data[f][W]['XSA1']*kq+data[f][W]['XSA2']*kq**2
    return R*data['kl'][W]['XSSM']- data['kl'][W]['XSSM']


LambdaNP=1e+3  #GeV
v= 246.
mh=125.1
mass ={
'ku':2.2e-3,
'kd':4.7e-3,
'ks':95e-3,
'kc':1.275,
}

# Helper function to confert from and to SMEFT
def kqtoCqH(kq,op):
    if op=='ku'or op=='kd'or op=='ks' or  op=='kc':
        return LambdaNP**2/v**3*(np.sqrt(2.0)*mass[op]*(1-kq))
    else:
        return kq

def CqHtokq(CqH,op):
    if op=='ku'or op=='kd'or op=='ks' or  op=='kc':
        return -(CqH/np.sqrt(2)/LambdaNP**2 *v**3/mass[op])+1
    else:
        return CqH
def CHtokl(CH):
    return 1-2.0*CH*v**4/mh**2/LambdaNP**2

def kltoCH(kl):
    return LambdaNP**2/v**4*mh**2*0.5*(1-kl)
