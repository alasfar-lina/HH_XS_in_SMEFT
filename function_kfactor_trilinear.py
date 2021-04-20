v=246.0
MT= 173.2
MH=125.1
K ={
'tri': 2.11645,
'box':1.84673,
'interference':1.9928
} # k_factors for each topology
def Kfac(chhh):
    """
    returns the K-factor for 14TeV HH production as a function of kappa_lambda
    """
    num =   (13.814233120011506 - 9.742976918098801*chhh + 2.1164397785287985*chhh**2)
    den = (7.480369598044151 - 4.885057884518588*chhh + 1.*chhh**2)
    return num/den


if __name__ == '__main__':
    print(Kfac(1.0))
