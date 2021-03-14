#!/usr/bin/python
from tkinter import *
from numpy import *


w_cnv = 750
h_cnv = 750

win_1 = Tk()
win_1.geometry("800x800")

canvas_1 = Canvas(win_1, width=w_cnv, height=h_cnv, bg="white")
canvas_1.pack(pady=20)

m_list = (1000,100,10,1)
(g,dag,hg,kg) = m_list
l_list = (1000,100,10,1)
(mm,cm,dm,m) = l_list
t_list = (1000,100,10,1)
(ms,cs,ds,s) = t_list

unit_m = 1*kg
unit_l = 1*m
unit_t = 1*ds

unit_v = unit_l/unit_t
unit_F = unit_m*unit_l/(unit_t**2)

grav = 9.82
rho_air = 1.204
rho_water = 1000.0
rho_steel = 8000.0
rho_rubber = 920.0
rho_medium = rho_air


obj = []
con = []

time = 0
dt = 0.05

timer = canvas_1.create_text(100,10,text="t = "+str(time),font="arial")

linecol = canvas_1.create_line( 0,0,1,0, dash = (2,2) )

def timestep():
    # Increment time
    global time
    time = time + dt
    canvas_1.itemconfigure(timer,text="t = "+str(time))

    # Zero forces
    for i in range(0,len(obj)):
        obj[i].F_other = [0,0]
        obj[i].F_R = [0,0]
    # Update constraint forces
    for i in range(0,len(con)):
        con[i].update()
    # Collision check
    for i in range(0,len(obj)):
        for j in range(i+1,len(obj)):
            Colcheck(obj[i],obj[j])
    # Take a time step
    for i in range(0,len(obj)):
        obj[i].integrate()
    # Loop
    win_1.after(10,timestep)

def Rotvec2(v,ang):
    return [ cos(ang)*v[0] - sin(ang)*v[1] , sin(ang)*v[0] + cos(ang)*v[1] ]

def Safediv(up,down):
    try:
        return up/down
    #except TypeError:
    #    return [ up[0]/down , up[1]/down ]
    except ZeroDivisionError:
        return sign(up)*1000000

def inters_p2l(point,line): # Intersection between normal line and main line, where normal line is at its closest point to object
    cooint = [0,0]
    # Distance
    #cooint[0] = ( (point.coo0[1]+point.coo0[0]/line.k)-line.cooA[1]+line.cooA[0]*line.k )/(line.k+1/line.k) # 2 line intersection
    cooint[0] = Safediv( ( point.coo0[1]+Safediv(point.coo0[0],line.k)-line.cooA[1]+line.cooA[0]*line.k ),( line.k+Safediv(1,line.k) ) ) # 2 line intersection
    cooint[0] = ( (cooint[0]>=line.xmin) and (cooint[0]<=line.xmax) )*cooint[0] + ( (cooint[0]<line.xmin) )*line.xmin + ( (cooint[0]>line.xmax) )*line.xmax
    cooint[1] = line.k*(cooint[0]-line.cooA[0])+line.cooA[1] # intersection x into fixline equation
    #canvas_1.coords(linecol,cooint[0],cooint[1],point.coo0[0],point.coo0[1])
    # Return value
    return cooint

def Colcheck(o1,o2):
    cooint = [0,0]
    distint = 0
    nvec = [0,0]
    if isinstance(o1, Ball):
        if isinstance(o2, Ball):
            dist_p2p = sqrt( (o2.coo0[0]-o1.coo0[0])**2 + (o2.coo0[1]-o1.coo0[1])**2 ) # distance between points
            if dist_p2p<0.5*(o1.d+o2.d):
                cosangle = (   (o2.coo0[0]-o1.coo0[0])/dist_p2p   )
                sinangle = (   (o2.coo0[1]-o1.coo0[1])/dist_p2p   )
                o1.F_R = add( o1.F_R , [ 100000000*cosangle*(dist_p2p-0.5*(o1.d+o2.d)) , 100000000*sinangle*(dist_p2p-0.5*(o1.d+o2.d)) ] )
                #o2.F_R = add( o2.F_R , [ -100000000*cosangle*(dist_p2p-0.5*(o1.d+o2.d)) , -100000000*sinangle*(dist_p2p-0.5*(o1.d+o2.d)) ] )
                o2.F_R = -o1.F_R
                o1.coo1 = o1.coo1 + 0.5*(o1.F_R/o1.m)*dt**2
                o2.coo1 = o2.coo1 + 0.5*(o2.F_R/o2.m)*dt**2

        #if obj[4] is o1:
        #if 1 == 1:
        if isinstance(o2, Fixline):
            # Coordinates of intersection point
            cooint = inters_p2l(o1,o2)
            #print(cooint)
            # Distance between object and intersection point
            distint = sqrt( (cooint[0]-o1.coo0[0])**2 + (cooint[1]-o1.coo0[1])**2 ) # distance between point and fixline
            #print(distint)
            # Normal vector
            nvec = divide( [ o1.coo0[0]-cooint[0] , o1.coo0[1]-cooint[1] ] , sqrt( (o1.coo0[0]-cooint[0])**2 + (o1.coo0[1]-cooint[1])**2 ) )
            if distint < 0.5*o1.d:
                o1.F_R = multiply( nvec , 500000000*(0.5*o1.d-distint) )
                o1.coo1 = o1.coo1 + 0.5*(o1.F_R/o1.m)*dt**2

class CnstSprDamp:
    def __init__(self,k,len0,damping,A,B):
        self.A = A
        self.B = B
        self.len0 = len0

        self.AB = sqrt(   (self.A.coo0[0]-self.B.coo0[0])**2   +   (self.A.coo0[1]-self.B.coo0[1])**2   )
        self.cosangle = (   (self.B.coo0[0]-self.A.coo0[0])/self.AB   )
        self.sinangle = (   (self.B.coo0[1]-self.A.coo0[1])/self.AB   )

        self.eps = (self.AB-self.len0)/self.len0

        self.vrel = sqrt( (self.A.v0[0] - self.B.v0[0])**2 + (self.A.v0[1] - self.B.v0[1])**2 )
        self.coodist = [abs(self.B.coo0[0] - self.A.coo0[0]), abs(self.B.coo0[1] - self.A.coo0[1])]
        self.line = canvas_1.create_line(self.A.coo0[0],self.A.coo0[1],self.B.coo0[0],self.B.coo0[1])

        self.damping = damping
        self.F_damping = [self.damping*self.cosangle*self.vrel, self.damping*self.sinangle*self.vrel]        
        
        self.len = [self.cosangle*(self.AB-self.len0), self.sinangle*(self.AB-self.len0)]
        self.k = k
        self.update()

    def update(self):
        self.ABm1 = self.AB
        self.AB = sqrt(   (self.A.coo0[0]-self.B.coo0[0])**2   +   (self.A.coo0[1]-self.B.coo0[1])**2   )
        self.deltaAB = self.AB - self.ABm1
        self.cosangle = (   (self.B.coo0[0]-self.A.coo0[0])/self.AB   )
        self.sinangle = (   (self.B.coo0[1]-self.A.coo0[1])/self.AB   )

        self.eps = (self.AB-self.len0)/self.len0

        self.vrel = sqrt( (self.A.v0[0] - self.B.v0[0])**2 + (self.A.v0[1] - self.B.v0[1])**2 )
        self.coodist = [abs(self.B.coo0[0] - self.A.coo0[0]), abs(self.B.coo0[1] - self.A.coo0[1])]
        self.len = [self.cosangle*(self.AB-self.len0), self.sinangle*(self.AB-self.len0)]

        self.F_damping = [self.damping*self.cosangle*self.deltaAB, self.damping*self.sinangle*self.deltaAB]      

        self.A.F_other = add( self.A.F_other , [(self.k*self.len[0] + self.F_damping[0]), (self.k*self.len[1] + self.F_damping[1])] )
        self.B.F_other = add( self.B.F_other , [-(self.k*self.len[0] + self.F_damping[0]), -(self.k*self.len[1] + self.F_damping[1])] )

        # Draw
        canvas_1.coords(self.line,self.A.coo0[0],self.A.coo0[1],self.B.coo0[0],self.B.coo0[1])
        if self.eps > 0.15:
            canvas_1.itemconfigure(self.line,fill="red")
        else:
            if self.eps < -0.15:
                canvas_1.itemconfigure(self.line,fill="blue",width=(1-self.eps)**6)
            else:
                canvas_1.itemconfigure(self.line,fill="black")
            
class Fixp:
    def __init__(self,canvas,x,y):
        self.coo0 = [x, y]
        self.d = 8
        self.v0 = [0,0]
        self.F_other = [0,0]
        self.ball = canvas.create_oval(x-self.d/2,y-self.d/2,x+self.d/2,y+self.d/2,fill="green")
    def integrate(self):
        return
        #print(1)
    
class Ball:
    def __init__(self,canvas,x,y,d,rho):
        # Indata
        self.coo1 = [x, y]
        self.coo0 = self.coo1
        self.coom1 = self.coo1
        self.d = d
        self.A = 0.25*pi*d**2
        self.V = 0.16666*pi*d**3
        self.rho = rho
        self.v0 = [0.0, 0.0]
        self.C_D = 0.47
        # Inertia
        self.m = self.rho*self.V
        self.I = 0.16667*m*d*d
        # Forces
        self.F_other = [0.0,0.0]
        self.F_other_abs = round(sqrt(self.F_other[0]**2 + self.F_other[1]**2)/10000)
        self.F_D = [0.0, 0.0]
        self.F_R = [0.0, 0.0]
        # Draw
        self.canvas = canvas
        self.ball = canvas.create_oval(x-d/2,y-d/2,x+d/2,y+d/2,fill="white")
        self.arrow_F = canvas.create_line(self.coo1[0],self.coo1[1],self.F_other[0],self.F_other[1],arrow=LAST,fill="grey")
        self.text_F = canvas_1.create_text(self.coo0[0],self.coo0[1],text=str(self.F_other_abs),font=("arial",8),fill="grey")
        self.arrow_g = canvas.create_line(self.coo1[0],self.coo1[1],self.coo1[0],self.coo1[1]+grav*5,arrow=LAST,fill="grey")
        self.offset_text_coo = -20
        self.text_coo = canvas_1.create_text(self.coo0[0],self.coo0[1]+self.offset_text_coo,text="("+str(round(self.coo0[0]))+", "+str(round(self.coo0[1]))+")",font=("arial",8))

    def colcheck(self):
        if self.coo1[1] + self.d/2 >= h_cnv:
            self.F_R = -100000000*((self.coo1[1] + self.d/2) - h_cnv)
            self.coo1[1] = self.coo1[1] + 0.5*(self.F_R/self.m)*dt**2 #0.5 här eller ej?
            canvas_1.itemconfigure(self.ball,fill="black")
        else:
            if self.coo1[0] + self.d/2 >= w_cnv:
                self.F_R = -100000000*((self.coo1[0] + self.d/2) - w_cnv)
                self.coo1[0] = self.coo1[0] + 0.5*(self.F_R/self.m)*dt**2 #0.5 här eller ej?
                canvas_1.itemconfigure(self.ball,fill="black")
            else:
                if self.coo1[0] - self.d/2 <= 0:
                    self.F_R = -100000000*((self.coo1[0] - self.d/2) - 0)
                    self.coo1[0] = self.coo1[0] + 0.5*(self.F_R/self.m)*dt**2 #0.5 här eller ej?
                    canvas_1.itemconfigure(self.ball,fill="black")
                else:
                    canvas_1.itemconfigure(self.ball,fill="white")

    def draw(self):
        self.canvas.coords(self.ball,self.coo1[0]-self.d/2,self.coo1[1]-self.d/2,self.coo1[0]+self.d/2,self.coo1[1]+self.d/2)

        self.canvas.coords(self.arrow_F,self.coo1[0],self.coo1[1],self.coo1[0]+self.F_other[0]/100000,self.coo1[1]+(self.F_other[1])/100000)
        self.canvas.coords(self.text_F,self.coo1[0]+self.F_other[0]/100000,self.coo1[1]+(self.F_other[1])/100000)

        self.canvas.itemconfigure(self.text_F,text=str(self.F_other_abs))
        self.canvas.coords(self.arrow_g,self.coo1[0],self.coo1[1],self.coo1[0],self.coo1[1]+grav*5)

        self.canvas.coords(self.text_coo,self.coo0[0],self.coo0[1]+self.offset_text_coo)
        self.canvas.itemconfigure(self.text_coo,text="("+str(round(self.coo0[0]))+", "+str(round(self.coo0[1]))+")")

    def integrate(self):
        #self.coom1 = subtract( self.coo0, multiply(self.v0,dt) )
        self.coom1 = self.coo0
        #print(self.coom1)
        self.coo0 = self.coo1
        #print(self.coo0)
        self.v0 = [0,0]#divide(subtract(self.coo0,self.coo1), dt)
        self.coo1 = [ 2*self.coo0[0] - self.coom1[0] + (self.F_other[0]/self.m)*dt**2 , 2*self.coo0[1] - self.coom1[1] + (self.F_other[1]/self.m + grav)*dt**2 ]
        #print(self.coo1)
        self.F_other_abs = round(sqrt(self.F_other[0]**2 + self.F_other[1]**2)/10000)

        self.colcheck()
        self.draw() # ändra så att draw är en allmän funktion som kallas efter den globala kollisionskollen

        #ändra så att den tar hänsyn till studsen
        self.v0 = [(self.coo1[0] - self.coo0[0])/dt, (self.coo1[1] - self.coo0[1])/dt]

class Fixline:
    def __init__(self,canvas,xA,yA,xB,yB):
        # Indata
        
        self.F_other = [0,0]
        self.cooA = [xA, yA]
        self.cooB = [xB, yB]
        self.xmax = max(self.cooA[0],self.cooB[0])
        self.xmin = min(self.cooA[0],self.cooB[0])
        self.AB = sqrt( (xB-xA)**2 + (yB-yA)**2 )
        self.ang = arccos( (xB-xA)/self.AB )
        self.coomidp = [ self.cooA[0] + 0.5*(self.cooB[0]-self.cooA[0]) , self.cooA[1] + 0.5*(self.cooB[1]-self.cooA[1]) ]
        self.unitnormalvec = Rotvec2( [1,0] , self.ang-pi/2 )
        self.normalvec = add( Rotvec2( [20,0] , self.ang-pi/2 ) , self.coomidp )
        self.nvec = [0,0]
        self.k = Safediv( (self.cooB[1]-self.cooA[1]) , (self.cooB[0]-self.cooA[0]) )
        self.cooint = [0,0]
        self.distint = 0
        self.text_distint = canvas_1.create_text(0,0,text="",font=("arial",8))

        # Forces
        self.F_other = [0.0,0.0]
        # Draw
        self.canvas = canvas
        self.lineAB = canvas.create_line(self.cooA[0],self.cooA[1],self.cooB[0],self.cooB[1])
        self.lineN = canvas.create_line(self.coomidp[0],self.coomidp[1],self.normalvec[0],self.normalvec[1],arrow=LAST)
        #self.linecol = canvas.create_line( 0,0,1,0, dash = (2,2) )


    def integrate(self):
        pass

    # def Distto(self,coop): # Intersection between normal line and main line, where normal line is at its closest point to object

    #     # Distance
    #     self.cooint[0] = ( coop[1]+coop[0]/self.k-self.cooA[1]+self.cooA[0]*self.k )/(self.k+1/self.k) # 2 line intersection
    #     self.cooint[0] = ( (self.cooint[0]>=self.xmin) and (self.cooint[0]<=self.xmax) )*self.cooint[0] + ( (self.cooint[0]<=self.xmin) )*self.xmin + ( (self.cooint[0]>=self.xmax) )*self.xmax
    #     self.cooint[1] = self.k*(self.cooint[0]-self.cooA[0])+self.cooA[1] # intersection x into fixline equation
    #     self.distint = sqrt( (self.cooint[0]-coop[0])**2 + (self.cooint[1]-coop[1])**2 ) # distance between point and fixline

    #     # Normal vector
    #     self.nvec = divide( [ coop[0]-self.cooint[0] , coop[1]-self.cooint[1] ] , sqrt( (coop[0]-self.cooint[0])**2 + (coop[1]-self.cooint[1])**2 ) )
        
    #     #print("cooint: "+str(self.cooint)+ "          objcoo: "+str(coop)+ "          dist: "+str(self.distint))

    #     # Draw
    #     self.canvas.coords(self.linecol,self.cooint[0],self.cooint[1],coop[0],coop[1])
    #     self.canvas.coords(self.text_distint, self.cooint[0] + 0.5*(coop[0] - self.cooint[0]) , self.cooint[1] + 0.5*(coop[1] - self.cooint[1]) )
    #     self.canvas.itemconfigure(self.text_distint,text=round(self.distint))

    #     #print(self.cooint)
        
    #     # Return value
    #     return self.distint
    
class Trace:
    def __init__(self,canvas,target,freq):
        self.canvas = canvas
        self.target = target
        self.coom1 = target.coom1
        self.freq = freq
        self.counter = 0

    def integrate(self):
        self.counter = self.counter + 1
        if self.counter >= self.freq:
            self.canvas.create_line( self.coom1[0],self.coom1[1],self.target.coo0[0],self.target.coo0[1],fill="grey", dash = (2,2) )
            self.coom1 = self.target.coom1
            self.counter = 0


obj.append(Ball(canvas_1,250,350,10,rho_rubber))
obj[0].v0[0] = 0.00
obj[0].v0[1] = -0.0

obj.append(Ball(canvas_1,300,250,10,rho_rubber))
obj[1].v0[0] = 0.0
obj[1].v0[1] = 0.0

obj.append(Ball(canvas_1,450,200,10,rho_rubber))
obj[2].v0[0] = 0.0
obj[2].v0[1] = 0.0

obj.append(Ball(canvas_1,470,30,10,rho_rubber))
obj[3].v0[0] = 0.0
obj[3].v0[1] = 0.0

obj.append(Ball(canvas_1,480,130,15,rho_rubber))
obj[4].v0[0] = 3.0
obj[4].v0[1] = 0.0

obj.append(Ball(canvas_1,130,50,15,rho_steel))
obj[5].v0[0] = 0.0
obj[5].v0[1] = 0.0

obj.append(Ball(canvas_1,100,15,8,rho_rubber))
obj[6].v0[0] = 0.0
obj[6].v0[1] = 0.0

obj.append(Ball(canvas_1,120,55,10,rho_rubber))
obj[7].v0[0] = 0.0
obj[7].v0[1] = 0.0

obj.append(Fixp(canvas_1,200,230))

obj.append(Ball(canvas_1,140,70,12,rho_rubber))
obj[9].v0[0] = 0.0
obj[9].v0[1] = 0.0

obj.append(Fixline(canvas_1,350,350,475,525))
obj.append(Fixline(canvas_1,450,500,600,550))
obj.append(Fixline(canvas_1,0,h_cnv-20,w_cnv,h_cnv-15)) # funkar ej om linjen är helt vågrät - FIXA

obj.append(Trace(canvas_1,obj[5],20))
obj.append(Fixp(canvas_1,400,400))
obj.append(Ball(canvas_1,400,450,30,rho_rubber))

con.append(CnstSprDamp(1000000,100,10000000,obj[0],obj[1]))
con.append(CnstSprDamp(1000000,100,10000000,obj[1],obj[2]))
con.append(CnstSprDamp(1000000,100,10000000,obj[2],obj[3]))
con.append(CnstSprDamp(1000000,100,10000000,obj[3],obj[0]))
con.append(CnstSprDamp(1000000,100,10000000,obj[2],obj[0]))
con.append(CnstSprDamp(1000000,100,10000000,obj[1],obj[3]))

con.append(CnstSprDamp(1000000,40,10000000,obj[5],obj[8]))

con.append(CnstSprDamp(1000000,25,10000000,obj[8],obj[6]))
con.append(CnstSprDamp(1000000,125,10000000,obj[6],obj[7]))
con.append(CnstSprDamp(1000000,125,10000000,obj[9],obj[6]))

con.append(CnstSprDamp(10000000,75,500000000,obj[14],obj[15]))

x_m = 0
y_m = 0
def motion(event):
    global x_m, y_m
    x_m, y_m = event.x, event.y
    #print('{}, {}'.format(x_m, y_m))
    obj[14].coo0 = [x_m, y_m]
    canvas_1.coords(obj[14].ball,x_m-obj[14].d/2,y_m-obj[14].d/2,x_m+obj[14].d/2,y_m+obj[14].d/2)
canvas_1.bind('<Motion>', motion)

#con.append(CnstSprDamp(1000000,100,10000000,obj[3],obj[0]))
#obj[5].coo0 = [x_m, y_m]
#obj[5].coo1 = [x_m, y_m]


## ATT GÖRA:
# fixa initial velocity för ball
# normalvectorn och enhetsnormalvectorn ballar ur för linjen, i vissa vinklar. beror nog på sinus periodicitet
# K fixa så att kontakt med linjen ej sker utanför dess kanter
# K ball to ball contact
# K contact with fixline
# friction for fixline OCH moment of inertia för ball
# cube 
# K skapa ett trace-objekt
# för fyrkanter/avlånga objekt osv, beräkna dess tvärsnittsarea som möter vinden varje tidssteg - påverkar luftmotsståndet - kan skapa vingar osv
# skapa musdrag av prylar
# fjäder/dämparna verkar ibland vara tjocka, trots att de är röda - fixa
# skapa rep, med inställning av antalet delar och längd
# allmän kodförbättring - förbättra integrationen så att den är mer samlad?
# återinför luftmotstånd OCH ljusblåa luftmotståndspilar
# skapa portaler?


timestep()

win_1.mainloop()
