import numpy as np
import matplotlib.pyplot as plt


r = np.random.rand(100000)

phir = 2*np.pi*r


g2 = 0.85
g= 0
theta = np.arange(0,181,1)
theta = theta * np.pi/180
mu = np.cos(theta)



PR = 3./4*(1+np.power(mu,2))
#PRmur = 3./4*(1+np.power(mur,2))

q = -8*r+4
D = 1+q*q/4
u = np.power((-q/2+np.sqrt(D)),1/3.)
mur = u-1/u




murHG = (np.power((2*r/(1-g2*g2)+1/(g2+1)),-2)-g2*g2-1)/(-2*g2)
PHGmur = 1-g*g/np.power((1+g*g-2*g*murHG),3/2.)
PHG = 1-g*g/np.power((1+g*g-2*g*mu),3/2.)


taur = -np.log(1-r)

thetar = np.arcsin(r)

#plt.subplot(2,1,1)
#plt.plot(mu,PR)

#plt.subplot(2,1,2)
#plt.hist(mur,20)

#for i in range(5):
    
    #g=-1 + 0.5*i
    
    #plt.subplot(5,1,i+1)
    #plt.plot(mu,1-g*g/np.power((1+g*g-2*g*mu),3/2.))
    

#plt.tight_layout()
#plt.show()

plt.subplot(4,1,1)
plt.hist(mur,20)
plt.subplot(4,1,2)
plt.hist(murHG,20)
plt.subplot(4,1,3)
plt.hist(taur,20)
plt.subplot(4,1,4)
plt.hist(thetar,20)

plt.tight_layout()
plt.show()


