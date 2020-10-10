
import numpy as np

def CLT(myE11,myE22,myG12,myNu12,myLaminate,myLaminateThickness,myPlyNumber):
    #--------------------------------------------------------------------------
    #
    # This function calculates the ABD Composite Stiffness Matrix Components
    # Classical Lamination Theory (CLT)
    #
    #--------------------------------------------------------------------------

    i=0    
    myThickness = myPlyNumber*myLaminateThickness
    Z = []
    Z.append(-myThickness/2.0)
    
    for j in range(1,myPlyNumber+1,1):
        Z.append(Z[j-1]+myLaminateThickness)
    
    # Calcuate Reduced Stiffnesses
    
    myQ11 = []
    myQ12 = []
    myQ16 = []
    myQ22 = []
    myQ26 = []
    myQ66 = []
    
    myQ11_S = []
    myQ12_S = []
    myQ16_S = []
    myQ22_S = []
    myQ26_S = []
    myQ66_S = []
    
    for i in range(0,1,1):
        myQ11.append((myE11**2)/(myE11-myNu12**2*myE22))
        myQ12.append((myNu12*myE11*myE22)/(myE11-myNu12**2*myE22))
        myQ16.append(0)
        myQ22.append((myE11*myE22)/(myE11-myNu12**2*myE22))
        myQ26.append(0)
        myQ66.append(myG12)
    
    i=0    
    for j in range(0,myPlyNumber,1):
        myQ11_S.append(myQ11[i]*np.cos(myLaminate[j]*np.pi/180.0)**4 + 2*(myQ12[i]+2*myQ66[i])*(np.cos(myLaminate[j]*np.pi/180.0))**2*(np.sin(myLaminate[j]*np.pi/180.0))**2 + myQ22[i]*(np.sin(myLaminate[j]*np.pi/180.0))**4)
        myQ12_S.append(myQ12[i]*((np.cos(myLaminate[j]*np.pi/180.0)**4)+(np.sin(myLaminate[j]*np.pi/180.0)**4))+ (myQ11[i]+myQ22[i]-4*myQ66[i])*np.cos(myLaminate[j]*np.pi/180.0)**2*np.sin(myLaminate[j]*np.pi/180.0)**2)
        myQ16_S.append((myQ11[i]-myQ12[i]-2*myQ66[i])*np.cos(myLaminate[j]*np.pi/180.0)**3*np.sin(myLaminate[j]*np.pi/180.0) - (myQ22[i] - myQ12[i] - 2*myQ66[i])*np.cos(myLaminate[j]*np.pi/180.0)*np.sin(myLaminate[j]*np.pi/180.0)**3)
        myQ22_S.append(myQ11[i]*np.sin(myLaminate[j]*np.pi/180.0)**4 + 2*(myQ12[i]+2*myQ66[i])*(np.cos(myLaminate[j]*np.pi/180.0))**2*(np.sin(myLaminate[j]*np.pi/180.0))**2 + myQ22[i]*(np.cos(myLaminate[j]*np.pi/180.0))**4)
        myQ26_S.append((myQ11[i]-myQ12[i]-2*myQ66[i])*np.cos(myLaminate[j]*np.pi/180.0)*np.sin(myLaminate[j]*np.pi/180.0)**3 - (myQ22[i] - myQ12[i] - 2*myQ66[i])*np.cos(myLaminate[j]*np.pi/180.0)**3*np.sin(myLaminate[j]*np.pi/180.0))
        myQ66_S.append((myQ11[i] + myQ22[i] - 2*myQ12[i] - 2*myQ66[i])*np.cos(myLaminate[j]*np.pi/180.0)**2*np.sin(myLaminate[j]*np.pi/180.0)**2 + myQ66[i]*(np.cos(myLaminate[j]*np.pi/180.0)**4+np.sin(myLaminate[j]*np.pi/180.0)**4))
    
    
    
    # Calcualte A Matrix
    #-------------
            
    A11_v = []
    A12_v = []
    A16_v = []
    A22_v = []
    A26_v = []
    A66_v = []
    i = 1
    for j in range(0,myPlyNumber,1):
        A11_v.append(myQ11_S[j]*(Z[i]-Z[i-1]))
        A12_v.append(myQ12_S[j]*(Z[i]-Z[i-1]))
        A16_v.append(myQ16_S[j]*(Z[i]-Z[i-1]))
        A22_v.append(myQ22_S[j]*(Z[i]-Z[i-1]))
        A26_v.append(myQ26_S[j]*(Z[i]-Z[i-1]))
        A66_v.append(myQ66_S[j]*(Z[i]-Z[i-1]))
        i = i+1
    
    A11 = sum(A11_v)
    A12 = sum(A12_v)
    A16 = sum(A16_v)
    A22 = sum(A22_v)
    A26 = sum(A26_v)
    A66 = sum(A66_v)
    
    # Calcualte B Matrix
    #-------------
            
    B11_v = []
    B12_v = []
    B16_v = []
    B22_v = []
    B26_v = []
    B66_v = []
    i = 1
    for j in range(0,myPlyNumber,1):
        B11_v.append(0.5*myQ11_S[j]*(Z[i]**2-Z[i-1]**2))
        B12_v.append(0.5*myQ12_S[j]*(Z[i]**2-Z[i-1]**2))
        B16_v.append(0.5*myQ16_S[j]*(Z[i]**2-Z[i-1]**2))
        B22_v.append(0.5*myQ22_S[j]*(Z[i]**2-Z[i-1]**2))
        B26_v.append(0.5*myQ26_S[j]*(Z[i]**2-Z[i-1]**2))
        B66_v.append(0.5*myQ66_S[j]*(Z[i]**2-Z[i-1]**2))
        i = i+1
    
    
    
    
    B11 = sum(B11_v)
    B12 = sum(B12_v)
    B16 = sum(B16_v)
    B22 = sum(B22_v)
    B26 = sum(B26_v)
    B66 = sum(B66_v)
    
    # Calcualte D Matrix
    #-------------
            
    D11_v = []
    D12_v = []
    D16_v = []
    D22_v = []
    D26_v = []
    D66_v = []
    i = 1
    for j in range(0,myPlyNumber,1):
        D11_v.append((1/3.0)*myQ11_S[j]*(Z[i]**3-Z[i-1]**3))
        D12_v.append((1/3.0)*myQ12_S[j]*(Z[i]**3-Z[i-1]**3))
        D16_v.append((1/3.0)*myQ16_S[j]*(Z[i]**3-Z[i-1]**3))
        D22_v.append((1/3.0)*myQ22_S[j]*(Z[i]**3-Z[i-1]**3))
        D26_v.append((1/3.0)*myQ26_S[j]*(Z[i]**3-Z[i-1]**3))
        D66_v.append((1/3.0)*myQ66_S[j]*(Z[i]**3-Z[i-1]**3))
        i = i+1
    
    
    D11 = sum(D11_v)
    D12 = sum(D12_v)
    D16 = sum(D16_v)
    D22 = sum(D22_v)
    D26 = sum(D26_v)
    D66 = sum(D66_v)
    return A11,A12,A16,A22,A26,A66,B11,B12,B16,B22,B26,B66,D11,D12,D16,D22,D26,D66


#-----------------------------------------------------------------------------------------------

# def CLT(myE11,myE22,myG12,myNu12,myLaminate,myLaminateThickness,myPlyNumber):






myE11 = 125744                       # Laminate E1 [MPa]  
myE22 = 10030                        # Laminate E2 [MPa]
myNu12 = 0.271                       # Laminate nu12 [-]   
myG12 = 5555                         # Laminate G12 [MPa]                                                                  
        


myLaminate = [45,-45,0,90,90,0,-45,45]
myPlyNumber = len(myLaminate)
myLaminateThickness = 0.125


A11,A12,A16,A22,A26,A66,B11,B12,B16,B22,B26,B66,D11,D12,D16,D22,D26,D66 = CLT(myE11,myE22,myG12,myNu12,myLaminate,myLaminateThickness,myPlyNumber)
        
        
