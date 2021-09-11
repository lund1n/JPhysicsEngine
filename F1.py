#!/usr/bin/python
from tkinter import *
from numpy import *
import time


w_cnv = 750
h_cnv = 750

win_1 = Tk()
win_1.geometry("800x800")

canvas_1 = Canvas(win_1, width=w_cnv, height=h_cnv, bg="white")
canvas_1.pack(pady=20)

#m_list = (1000,100,10,1)
#(g,dag,hg,kg) = m_list
#l_list = (1000,100,10,1)
#(mm,cm,dm,m) = l_list
#t_list = (1000,100,10,1)
#(ms,cs,ds,s) = t_list

#unit_m = 1*kg
#unit_l = 1*m
#unit_t = 1*ds

#unit_v = unit_l/unit_t
#unit_F = unit_m*unit_l/(unit_t**2)

grav = [ 0 , 9.82 ]
rho_air = 1.204
rho_water = 1000.0
rho_steel = 8000.0
rho_rubber = 920.0
rho_medium = rho_air

obj = []
con = []
viz_sf = 40
viz_sf2 = 1000000
sf_farrow1 = 200
sf_farrow2 = 100000000
viz = []

#colcounter = 0

t_elapsed = 0
dt = 0.025

dtr = canvas_1.create_text(10,20,text="dt = "+str(dt),font="arial",anchor=SW)
t_elapsedr = canvas_1.create_text(10,35,text="t = "+str(round(t_elapsed,9)),font="arial",anchor=SW)

linecol = canvas_1.create_line( 0,0,1,0, dash = (2,2) )

def Timestep():
    # Increment t_elapsed
    global t_elapsed
    t_elapsed = t_elapsed + dt
    canvas_1.itemconfigure(t_elapsedr,text="t = "+str(round(t_elapsed,9)))
    #time.sleep(0.2)

    # Remove collision visualisation
    for i in range(0,len(viz)):
        canvas_1.delete(viz[0])
        viz.pop(0)
    # Zero forces
    for i in range(0,len(obj)):
        obj[i].F_other = [0,0]
        obj[i].tau_other = 0.0
        obj[i].F_R = [0,0]
    # Update constraint forces
    for i in range(0,len(con)):
        con[i].update()
    # Collision check
    for i in range(0,len(obj)):
        for j in range(i+1,len(obj)):
            Check_CollisionType(obj[i],obj[j])
    # Take a time step
    for i in range(0,len(obj)):
        obj[i].integrate()
        #if i==4:
            #print(obj[i].theta1)
            #print( "Vel:   "+str((obj[i].coo1[1]-obj[i].coo0[1])/dt) )
            #print((obj[i].coo1[1]-obj[i].coo0[1])/dt)
            # if (sign(obj[i].coo1[1]-obj[i].coo0[1])) != (sign(obj[i].coo0[1]-obj[i].coom1[1])):
            #     canvas_1.create_oval(obj[i].coo1[0]-obj[i].r,obj[i].coo1[1]-obj[i].r,obj[i].coo1[0]+obj[i].r,obj[i].coo1[1]+obj[i].r,outline="blue")
            #     print(obj[i].coo1)
        #if i == 4:
            #print("Object_Ball:     " + str(obj[i].coo1))
            #print( sqrt( (obj[4].coo0[0]-obj[4].coom1[0])**2+(obj[4].coo0[1]-obj[4].coom1[1])**2 ) )
        #if i == 36:
            #print(obj[i].F_D)
        #if i == 15: # the mouse ball
        #    print(obj[i].F_D)
    
    #print(obj[4].coo0[1]-obj[4].coom1[1])
    #print( sqrt( (obj[4].coo0[0]-obj[4].coom1[0])**2+(obj[4].coo0[1]-obj[4].coom1[1])**2 ) )

    # Loop
    win_1.after(10,Timestep)

def Rotvec2(v,ang):
    return array([ cos(ang)*v[0] - sin(ang)*v[1] , sin(ang)*v[0] + cos(ang)*v[1] ])

def Rotvec2_90degcw(v):
    return array([ v[1] , -v[0] ])

def Rotzvec3_90degcw(v):
    return array([ v[1] , -v[0] , v[2] ])

def Safediv(up,down):
    if down != 0:
        return up/down
    #except TypeError:
    #    return [ up[0]/down , up[1]/down ]
    else:
        return sign(up)*100000000

def Crossprod_vec2(v,w):
    return v[0]*w[1]-v[1]*w[0]

def Crossprod_vec3(v,w):
    return [ v[1]*w[2] - v[2]*w[1] , v[2]*w[0] - v[0]*w[2] , v[0]*w[1] - v[1]*w[0] ]

def Matmultcol(M,v):
    return [M[0][0]*v[0] , M[1][1]*v[1] , M[2][2]*v[2]]

def Scalarproj_vec3(v,w):
    return (v[0]*w[0] + v[1]*w[1] + v[2]*w[2])/sqrt(w[0]**2 + w[1]**2 + w[2]**2)

def Vecproj_vec2(v,w):
    return multiply( Safediv( v[0]*w[0] + v[1]*w[1] , (w[0]**2 + w[1]**2) ) , w )

def Vecproj_vec3(v,w):
    return multiply( Safediv( v[0]*w[0] + v[1]*w[1] + v[2]*w[2] , (w[0]**2 + w[1]**2 + w[2]**2) ) , w )

def ClosestPointOnLineSegment(px,py,line):
    cooint = [0,0]
    # Distance
    # Intersection between normal line and main line, where normal line is at its closest point to object
    # 2 line perpendicular intersection (with exception for line k-value = 0 <-> horizontal line )
    cooint[0] = (line.k!=0)*Safediv( ( py+Safediv(px,line.k)-line.cooA[1]+line.cooA[0]*line.k ),( line.k+Safediv(1,line.k) ) ) + (line.k==0)*px
    # If the point is outside of the line boundaries, set the x coordinate to the corresponding line endpoint
    cooint[0] = ( (cooint[0]>=line.coomin[0]) and (cooint[0]<=line.coomax[0]) )*cooint[0] + ( (cooint[0]<line.coomin[0]) )*line.coomin[0] + ( (cooint[0]>line.coomax[0]) )*line.coomax[0]
    # Set the y value and SPECIAL CASE: (with exception for line endpoint x coordinates being the same <-> vertical line & in that case, check whether the point is inbetween the line endpoint y coordinates)
    cooint[1] = (line.cooA[0]!=line.cooB[0])*( line.k*(cooint[0]-line.cooA[0])+line.cooA[1] )   +   (line.cooA[0]==line.cooB[0])*(   (py<=line.coomax[1] and py>=line.coomin[1])*py   +   (py>line.coomax[1])*line.coomax[1]   +   (py<line.coomin[1])*line.coomin[1]   ) # intersection x into fixline equation
    #canvas_1.coords(linecol,cooint[0],cooint[1],point.coo0[0],point.coo0[1])
    # Return value
    return cooint

def ClosestPointOnLineSegmentEdgeIndicator(px,py,line):
    cooint = [0,0,0]
    # Distance
    # Intersection between normal line and main line, where normal line is at its closest point to object
    # 2 line perpendicular intersection (with exception for line k-value = 0 <-> horizontal line )
    cooint[0] = (line.k!=0)*Safediv( ( py+Safediv(px,line.k)-line.cooA[1]+line.cooA[0]*line.k ),( line.k+Safediv(1,line.k) ) ) + (line.k==0)*px
    cooint[2] = (cooint[0]<line.coomin[0]) or (cooint[0]>line.coomax[0]) or (line.cooA[0]==line.cooB[0])*( (py>line.coomax[1]) or (py<line.coomin[1]) )
    # If the point is outside of the line boundaries, set the x coordinate to the corresponding line endpoint
    cooint[0] = ( (cooint[0]>=line.coomin[0]) and (cooint[0]<=line.coomax[0]) )*cooint[0] + ( (cooint[0]<line.coomin[0]) )*line.coomin[0] + ( (cooint[0]>line.coomax[0]) )*line.coomax[0]
    # Set the y value and SPECIAL CASE: (with exception for line endpoint x coordinates being the same <-> vertical line & in that case, check whether the point is inbetween the line endpoint y coordinates)
    cooint[1] = (line.cooA[0]!=line.cooB[0])*( line.k*(cooint[0]-line.cooA[0])+line.cooA[1] )   +   (line.cooA[0]==line.cooB[0])*(   (py<=line.coomax[1] and py>=line.coomin[1])*py   +   (py>line.coomax[1])*line.coomax[1]   +   (py<line.coomin[1])*line.coomin[1]   ) # intersection x into fixline equation
    #canvas_1.coords(linecol,cooint[0],cooint[1],point.coo0[0],point.coo0[1])
    # Return value
    return cooint

def IntersectionPoint_Line_Line(l1,l2):
    cooint = [0,0]
    # 2 line intersection
    cooint[0] = Safediv( ( (l2.cooA[1]-l2.k*l2.cooA[0])-(l1.cooA[1]-l1.k*l1.cooA[0]) ),(l1.k-l2.k) ) #xi = (m2-m1)/(k1-k2) = ( (y2-k2x2)-(y1-k1x1) )/(k1-k2)
    cooint[1] = l1.k*cooint[0]+(l1.cooA[1]-l1.k*l1.cooA[0]) # yi = k1xi+m1 = k1xi+(y1-k1x1)
    return cooint

def Check_CollisionType(o1,o2):
    if isinstance(o1, Object_Ball):
        if isinstance(o2, Object_Ball):
            Collision_Ball_Ball(o1,o2)
            return
        if isinstance(o2, Object_Line):
            Collision_Ball_Line(o1,o2)
        #if obj[4] is o1:
        #if 1 == 1:
        if isinstance(o2, Object_FixedLine):
            Collision_Ball_FixedLine(o1,o2)
            return
    if isinstance(o1, Object_FixedLine):
        if isinstance(o2, Object_Ball):
            Collision_Ball_FixedLine(o2,o1)
            return
        if isinstance(o2, Object_Line):
            Collision_Line_FixedLine(o2,o1)
            return
    if isinstance(o1, Object_Line):
        if isinstance(o2, Object_FixedLine):
            Collision_Line_FixedLine(o1,o2)
        if isinstance(o2, Object_Ball):
            Collision_Ball_Line(o2,o1)
        if isinstance(o2, Object_Line):
            Collision_Line_Line(o2,o1)
            return

def Contact_line_line(o1,o2,coocol,unitvec_n,unitvec_t):

    rcol_b1 = [coocol[0]-o1.coo1[0],coocol[1]-o1.coo1[1],0]
    rcol_b2 = [coocol[0]-o2.coo1[0],coocol[1]-o2.coo1[1],0]
    vel_o1 = divide( [ o1.coo1[0]-o1.coo0[0] , o1.coo1[1]-o1.coo0[1] , 0 ] , dt )
    vel_o2 = divide( [ o2.coo1[0]-o2.coo0[0] , o2.coo1[1]-o2.coo0[1] , 0 ] , dt )
    omega0_o1 = (o1.theta1-o1.theta0)/dt
    omega0_o2 = (o2.theta1-o2.theta0)/dt

    cpo1 = Crossprod_vec3( [0,0,omega0_o1] , [rcol_b1[0],rcol_b1[1],0] )
    cpo2 = Crossprod_vec3( [0,0,omega0_o2] , [rcol_b2[0],rcol_b2[1],0] )
    vel0_rel = subtract( add( vel_o2 , [cpo2[0],cpo2[1],0] ) , add( vel_o1 , [cpo1[0],cpo1[1],0] ) )

    upper_n = -(1+o1.C_el*o2.C_el)*Scalarproj_vec3( vel0_rel,unitvec_n)
    #lower_n = (o1.invm + o2.invm + Scalarproj_vec3(Matmultcol( o1.invImat , Crossprod_vec3(Crossprod_vec3(rcol_o1,unitvec_n),rcol_o1) ) + Matmultcol( o2.invImat , Crossprod_vec3(Crossprod_vec3(rcol_o2,unitvec_n),rcol_o2) ) ,unitvec_n) )
    lower_n = (o1.invm + o2.invm + divide( (Crossprod_vec2(rcol_b1,unitvec_n))**2 , o1.I ) + divide( (Crossprod_vec2(rcol_b2,unitvec_n))**2 , o2.I ) )


    F_n = (upper_n/lower_n)/dt

    o1.F_other =  subtract( o1.F_other , multiply([unitvec_n[0],unitvec_n[1]],F_n) )
    o2.F_other =  add( o2.F_other , multiply([unitvec_n[0],unitvec_n[1]],F_n) )

    upper_t = -(o1.C_fric*o2.C_fric)*Scalarproj_vec3( vel0_rel,unitvec_t)
    #lower_t = (o1.invm + o2.invm + Scalarproj_vec3(Matmultcol( o1.invImat , Crossprod_vec3(Crossprod_vec3(rcol_o1,unitvec_t),rcol_o1) ) + Matmultcol( o2.invImat , Crossprod_vec3(Crossprod_vec3(rcol_o2,unitvec_t),rcol_o2) ) ,unitvec_t) )
    lower_t = (o1.invm + o2.invm + divide( (Crossprod_vec2(rcol_b1,unitvec_t))**2 , o1.I ) + divide( (Crossprod_vec2(rcol_b2,unitvec_t))**2 , o2.I ) )


    F_t = (upper_t/lower_t)/dt

    tauo1_n = multiply(F_n,Crossprod_vec3(rcol_b1,unitvec_n))
    tauo2_n = multiply(F_n,Crossprod_vec3(rcol_b2,unitvec_n))
    o1.tau_other =  subtract( o1.tau_other , tauo1_n[2] )
    o2.tau_other =  add( o2.tau_other , tauo2_n[2] )

    tauo1_t = multiply(F_t,Crossprod_vec3(rcol_b1,unitvec_t))
    tauo2_t = multiply(F_t,Crossprod_vec3(rcol_b2,unitvec_t))
    o1.tau_other =  subtract( o1.tau_other , tauo1_t[2] )
    o2.tau_other =  add( o2.tau_other , tauo2_t[2] )

    o1.F_other =  subtract( o1.F_other , multiply([unitvec_t[0],unitvec_t[1]],F_t) )
    o2.F_other =  add( o2.F_other , multiply([unitvec_t[0],unitvec_t[1]],F_t) )


    # Draw contact point
    viz.append( canvas_1.create_oval(coocol[0]-4,coocol[1]-4,coocol[0]+4,coocol[1]+4,fill="red",outline="red") )
    # Draw force arrows
    #viz.append( canvas_1.create_line(coocol[0]-unitvec_n[0]*viz_sf,coocol[1]-unitvec_n[1]*viz_sf,coocol[0]+unitvec_n[0]*viz_sf,coocol[1]+unitvec_n[1]*viz_sf,arrow=BOTH,fill="green") )
    viz.append( canvas_1.create_line(coocol[0]-unitvec_n[0]*sf_farrow1*tanh(abs(F_n)/sf_farrow2),coocol[1]-unitvec_n[1]*sf_farrow1*tanh(abs(F_n)/sf_farrow2),coocol[0]+unitvec_n[0]*sf_farrow1*tanh(abs(F_n)/sf_farrow2),coocol[1]+unitvec_n[1]*sf_farrow1*tanh(abs(F_n)/sf_farrow2),arrow=BOTH,fill="red") )
    # Draw tang vel arrows
    #viz.append( canvas_1.create_line(coocol[0]-unitvec_t[0]*30,coocol[1]-unitvec_t[1]*30,coocol[0]+unitvec_t[0]*30,coocol[1]+unitvec_t[1]*30,arrow=BOTH,fill="green") )
    viz.append( canvas_1.create_line(coocol[0]-unitvec_t[0]*sf_farrow1*tanh(abs(F_t)/sf_farrow2),coocol[1]-unitvec_t[1]*sf_farrow1*tanh(abs(F_t)/sf_farrow2),coocol[0]+unitvec_t[0]*sf_farrow1*tanh(abs(F_t)/sf_farrow2),coocol[1]+unitvec_t[1]*sf_farrow1*tanh(abs(F_t)/sf_farrow2),arrow=BOTH,fill="red") )

def Contact_dyn_statline(o1,o2,coocol,unitvec_n,unitvec_t):
    rcol_o1 = [coocol[0]-o1.coo1[0],coocol[1]-o1.coo1[1],0]

    vel_o1 = divide( [ o1.coo1[0]-o1.coo0[0] , o1.coo1[1]-o1.coo0[1] , 0 ] , dt )
    omega0_o1 = (o1.theta1-o1.theta0)/dt
    cpo1 = Crossprod_vec3( [0,0,omega0_o1] , [rcol_o1[0],rcol_o1[1],0] )
    vel0_rel = add( vel_o1 , [cpo1[0],cpo1[1],0] )

    upper_n = -(1+o1.C_el*o2.C_el)*Scalarproj_vec3( vel0_rel,unitvec_n)
    lower_n = (o1.invm + divide( (Crossprod_vec2(rcol_o1,unitvec_n))**2 , o1.I ) )
    F_n = -(upper_n/lower_n)/dt
    F_n_vec = multiply([unitvec_n[0],unitvec_n[1]],F_n)

    upper_t = -(o1.C_fric*o2.C_fric)*Scalarproj_vec3( vel0_rel,unitvec_t)
    lower_t = (o1.invm + divide( (Crossprod_vec2(rcol_o1,unitvec_t))**2 , o1.I ) )
    F_t = -(upper_t/lower_t)/dt
    F_t_vec = multiply([unitvec_t[0],unitvec_t[1]],F_t)
    F_t_unitvec = Safediv(F_t_vec,-abs(F_t)) # Nödvändig pga att unitvec_t ej pekar åt rätt håll

    o1.F_other =  subtract( o1.F_other , F_n_vec )
    o1.F_other =  subtract( o1.F_other , F_t_vec )

    tauo1_n = multiply(F_n,Crossprod_vec3(rcol_o1,unitvec_n))
    tauo1_t = multiply(F_t,Crossprod_vec3(rcol_o1,unitvec_t))

    o1.tau_other =  subtract( o1.tau_other , tauo1_n[2] )
    o1.tau_other =  subtract( o1.tau_other , tauo1_t[2] )
    
    # Draw contact point
    viz.append( canvas_1.create_oval(coocol[0]-4,coocol[1]-4,coocol[0]+4,coocol[1]+4,fill="red",outline="red") )
    # Draw normal force arrow
    #viz.append( canvas_1.create_line(coocol[0],coocol[1],coocol[0]+unitvec_n[0]*viz_sf,coocol[1]+unitvec_n[1]*viz_sf,arrow=LAST,fill="red") )
    #viz.append( canvas_1.create_line(coocol[0],coocol[1],coocol[0]+unitvec_n[0]*abs(F_n)/viz_sf2,coocol[1]+unitvec_n[1]*abs(F_n)/viz_sf2,arrow=LAST,fill="red") )
    viz.append( canvas_1.create_line(coocol[0],coocol[1],coocol[0]+unitvec_n[0]*sf_farrow1*tanh(abs(F_n)/sf_farrow2),coocol[1]+unitvec_n[1]*sf_farrow1*tanh(abs(F_n)/sf_farrow2),arrow=LAST,fill="red") )
    # Draw normal force arrow
    #viz.append( canvas_1.create_line(coocol[0],coocol[1],coocol[0]+F_t_unitvec[0]*viz_sf,coocol[1]+F_t_unitvec[1]*viz_sf,arrow=LAST,fill="red") )
    #viz.append( canvas_1.create_line(coocol[0],coocol[1],coocol[0]+F_t_unitvec[0]*abs(F_t)/viz_sf2,coocol[1]+F_t_unitvec[1]*abs(F_t)/viz_sf2,arrow=LAST,fill="red") )
    viz.append( canvas_1.create_line(coocol[0],coocol[1],coocol[0]+F_t_unitvec[0]*sf_farrow1*tanh(abs(F_t)/sf_farrow2),coocol[1]+F_t_unitvec[1]*sf_farrow1*tanh(abs(F_t)/sf_farrow2),arrow=LAST,fill="red") )

def Contact_dyn_line(o1,o2,coocol,unitvec_n,unitvec_t):
    
    rcol_o1 = [coocol[0]-o1.coo1[0],coocol[1]-o1.coo1[1],0]
    rcol_o2 = [coocol[0]-o2.coo1[0],coocol[1]-o2.coo1[1],0]
    vel_o1 = divide( [ o1.coo1[0]-o1.coo0[0] , o1.coo1[1]-o1.coo0[1] , 0 ] , dt )
    vel_o2 = divide( [ o2.coo1[0]-o2.coo0[0] , o2.coo1[1]-o2.coo0[1] , 0 ] , dt )
    omega0_o1 = (o1.theta1-o1.theta0)/dt
    omega0_o2 = (o2.theta1-o2.theta0)/dt
    #unitvec_n = [o2.coo1[0]-o1.coo1[0],o2.coo1[1]-o1.coo1[1],0]/dist_o1o2 # Unit normal vector
    #unitvec_t = Rotzvec3_90degcw(unitvec_n) # Unit tangential vector

    cpo1 = Crossprod_vec3( [0,0,omega0_o1] , [rcol_o1[0],rcol_o1[1],0] )
    cpo2 = Crossprod_vec3( [0,0,omega0_o2] , [rcol_o2[0],rcol_o2[1],0] )
    vel0_rel = subtract( add( vel_o2 , [cpo2[0],cpo2[1],0] ) , add( vel_o1 , [cpo1[0],cpo1[1],0] ) )

    upper_n = -(1+o1.C_el*o2.C_el)*Scalarproj_vec3( vel0_rel,unitvec_n)
    #lower_n = (o1.invm + o2.invm + Scalarproj_vec3(Matmultcol( o1.invImat , Crossprod_vec3(Crossprod_vec3(rcol_o1,unitvec_n),rcol_o1) ) + Matmultcol( o2.invImat , Crossprod_vec3(Crossprod_vec3(rcol_o2,unitvec_n),rcol_o2) ) ,unitvec_n) )
    lower_n = (o1.invm + o2.invm + divide( (Crossprod_vec2(rcol_o1,unitvec_n))**2 , o1.I ) + divide( (Crossprod_vec2(rcol_o2,unitvec_n))**2 , o2.I ) )


    F_n = (upper_n/lower_n)/dt

    o1.F_other =  subtract( o1.F_other , multiply([unitvec_n[0],unitvec_n[1]],F_n) )
    o2.F_other =  add( o2.F_other , multiply([unitvec_n[0],unitvec_n[1]],F_n) )

    upper_t = -(o1.C_fric*o2.C_fric)*Scalarproj_vec3( vel0_rel,unitvec_t)
    #lower_t = (o1.invm + o2.invm + Scalarproj_vec3(Matmultcol( o1.invImat , Crossprod_vec3(Crossprod_vec3(rcol_o1,unitvec_t),rcol_o1) ) + Matmultcol( o2.invImat , Crossprod_vec3(Crossprod_vec3(rcol_o2,unitvec_t),rcol_o2) ) ,unitvec_t) )
    lower_t = (o1.invm + o2.invm + divide( (Crossprod_vec2(rcol_o1,unitvec_t))**2 , o1.I ) + divide( (Crossprod_vec2(rcol_o2,unitvec_t))**2 , o2.I ) )


    F_t = (upper_t/lower_t)/dt

    tauo1_t = multiply(F_t,Crossprod_vec3(rcol_o1,unitvec_t))
    tauo2_t = multiply(F_t,Crossprod_vec3(rcol_o2,unitvec_t))
    o1.tau_other =  subtract( o1.tau_other , tauo1_t[2] )
    o2.tau_other =  add( o2.tau_other , tauo2_t[2] )

    tauo2_n = multiply(F_n,Crossprod_vec3(rcol_o2,unitvec_n))
    o2.tau_other =  add( o2.tau_other , tauo2_n[2] )

    o1.F_other =  subtract( o1.F_other , multiply([unitvec_t[0],unitvec_t[1]],F_t) )
    o2.F_other =  add( o2.F_other , multiply([unitvec_t[0],unitvec_t[1]],F_t) )
    
    # Draw contact point
    viz.append( canvas_1.create_oval(coocol[0]-4,coocol[1]-4,coocol[0]+4,coocol[1]+4,fill="red",outline="red") )
    # Draw force arrows
    #viz.append( canvas_1.create_line(coocol[0]-unitvec_n[0]*viz_sf,coocol[1]-unitvec_n[1]*viz_sf,coocol[0]+unitvec_n[0]*viz_sf,coocol[1]+unitvec_n[1]*viz_sf,arrow=BOTH,fill="green") )
    viz.append( canvas_1.create_line(coocol[0]-unitvec_n[0]*sf_farrow1*tanh(abs(F_n)/sf_farrow2),coocol[1]-unitvec_n[1]*sf_farrow1*tanh(abs(F_n)/sf_farrow2),coocol[0]+unitvec_n[0]*sf_farrow1*tanh(abs(F_n)/sf_farrow2),coocol[1]+unitvec_n[1]*sf_farrow1*tanh(abs(F_n)/sf_farrow2),arrow=BOTH,fill="red") )
    # Draw tang vel arrows
    #viz.append( canvas_1.create_line(coocol[0]-unitvec_t[0]*30,coocol[1]-unitvec_t[1]*30,coocol[0]+unitvec_t[0]*30,coocol[1]+unitvec_t[1]*30,arrow=BOTH,fill="green") )
    viz.append( canvas_1.create_line(coocol[0]-unitvec_t[0]*sf_farrow1*tanh(abs(F_t)/sf_farrow2),coocol[1]-unitvec_t[1]*sf_farrow1*tanh(abs(F_t)/sf_farrow2),coocol[0]+unitvec_t[0]*sf_farrow1*tanh(abs(F_t)/sf_farrow2),coocol[1]+unitvec_t[1]*sf_farrow1*tanh(abs(F_t)/sf_farrow2),arrow=BOTH,fill="red") )

def Collision_Ball_Ball(b1,b2):
    dist_b1b2 = sqrt( (b2.coo1[0]-b1.coo1[0])**2 + (b2.coo1[1]-b1.coo1[1])**2 ) # distance between points
    if dist_b1b2<(b1.r+b2.r):
        
        coocol = add(b1.coo1,(b1.r/dist_b1b2)*subtract(b2.coo1,b1.coo1)) # coordinates of collision point
        unitvec_n = [b2.coo1[0]-b1.coo1[0],b2.coo1[1]-b1.coo1[1],0]/dist_b1b2 # Unit normal vector
        unitvec_t = Rotzvec3_90degcw(unitvec_n) # Unit tangential vector
        
        rcol_b1 = [coocol[0]-b1.coo1[0],coocol[1]-b1.coo1[1],0]
        rcol_b2 = [coocol[0]-b2.coo1[0],coocol[1]-b2.coo1[1],0]
        vel_o1 = divide( [ b1.coo1[0]-b1.coo0[0] , b1.coo1[1]-b1.coo0[1] , 0 ] , dt )
        vel_o2 = divide( [ b2.coo1[0]-b2.coo0[0] , b2.coo1[1]-b2.coo0[1] , 0 ] , dt )
        omega0_o1 = (b1.theta1-b1.theta0)/dt
        omega0_o2 = (b2.theta1-b2.theta0)/dt

        cpo1 = Crossprod_vec3( [0,0,omega0_o1] , [rcol_b1[0],rcol_b1[1],0] )
        cpo2 = Crossprod_vec3( [0,0,omega0_o2] , [rcol_b2[0],rcol_b2[1],0] )
        vel0_rel = subtract( add( vel_o2 , [cpo2[0],cpo2[1],0] ) , add( vel_o1 , [cpo1[0],cpo1[1],0] ) )

        upper_n = -(1+b1.C_el*b2.C_el)*Scalarproj_vec3( vel0_rel,unitvec_n)
        #lower_n = (b1.invm + b2.invm + Scalarproj_vec3(Matmultcol( b1.invImat , Crossprod_vec3(Crossprod_vec3(rcol_o1,unitvec_n),rcol_o1) ) + Matmultcol( b2.invImat , Crossprod_vec3(Crossprod_vec3(rcol_o2,unitvec_n),rcol_o2) ) ,unitvec_n) )
        lower_n = (b1.invm + b2.invm + divide( (Crossprod_vec2(rcol_b1,unitvec_n))**2 , b1.I ) + divide( (Crossprod_vec2(rcol_b2,unitvec_n))**2 , b2.I ) )


        F_n = (upper_n/lower_n)/dt

        b1.F_other =  subtract( b1.F_other , multiply([unitvec_n[0],unitvec_n[1]],F_n) )
        b2.F_other =  add( b2.F_other , multiply([unitvec_n[0],unitvec_n[1]],F_n) )

        upper_t = -(b1.C_fric*b2.C_fric)*Scalarproj_vec3( vel0_rel,unitvec_t)
        #lower_t = (b1.invm + b2.invm + Scalarproj_vec3(Matmultcol( b1.invImat , Crossprod_vec3(Crossprod_vec3(rcol_o1,unitvec_t),rcol_o1) ) + Matmultcol( b2.invImat , Crossprod_vec3(Crossprod_vec3(rcol_o2,unitvec_t),rcol_o2) ) ,unitvec_t) )
        lower_t = (b1.invm + b2.invm + divide( (Crossprod_vec2(rcol_b1,unitvec_t))**2 , b1.I ) + divide( (Crossprod_vec2(rcol_b2,unitvec_t))**2 , b2.I ) )


        F_t = (upper_t/lower_t)/dt

        tauo1 = multiply(F_t,Crossprod_vec3(rcol_b1,unitvec_t))
        tauo2 = multiply(F_t,Crossprod_vec3(rcol_b2,unitvec_t))


        b1.tau_other =  subtract( b1.tau_other , tauo1[2] )
        b2.tau_other =  add( b2.tau_other , tauo2[2] )

        b1.F_other =  subtract( b1.F_other , multiply([unitvec_t[0],unitvec_t[1]],F_t) )
        b2.F_other =  add( b2.F_other , multiply([unitvec_t[0],unitvec_t[1]],F_t) )
    
        # Draw contact point
        viz.append( canvas_1.create_oval(coocol[0]-4,coocol[1]-4,coocol[0]+4,coocol[1]+4,fill="red",outline="red") )
        # Draw force arrows
        #viz.append( canvas_1.create_line(coocol[0]-unitvec_n[0]*viz_sf,coocol[1]-unitvec_n[1]*viz_sf,coocol[0]+unitvec_n[0]*viz_sf,coocol[1]+unitvec_n[1]*viz_sf,arrow=BOTH,fill="green") )
        viz.append( canvas_1.create_line(coocol[0]-unitvec_n[0]*sf_farrow1*tanh(abs(F_n)/sf_farrow2),coocol[1]-unitvec_n[1]*sf_farrow1*tanh(abs(F_n)/sf_farrow2),coocol[0]+unitvec_n[0]*sf_farrow1*tanh(abs(F_n)/sf_farrow2),coocol[1]+unitvec_n[1]*sf_farrow1*tanh(abs(F_n)/sf_farrow2),arrow=BOTH,fill="red") )
        # Draw tang vel arrows
        #viz.append( canvas_1.create_line(coocol[0]-unitvec_t[0]*30,coocol[1]-unitvec_t[1]*30,coocol[0]+unitvec_t[0]*30,coocol[1]+unitvec_t[1]*30,arrow=BOTH,fill="green") )
        viz.append( canvas_1.create_line(coocol[0]-unitvec_t[0]*sf_farrow1*tanh(abs(F_t)/sf_farrow2),coocol[1]-unitvec_t[1]*sf_farrow1*tanh(abs(F_t)/sf_farrow2),coocol[0]+unitvec_t[0]*sf_farrow1*tanh(abs(F_t)/sf_farrow2),coocol[1]+unitvec_t[1]*sf_farrow1*tanh(abs(F_t)/sf_farrow2),arrow=BOTH,fill="red") )


        # Move balls apart, without adding velocity
        pendist = (dist_b1b2-(b1.r+b2.r))
        b1.coo1 = add( b1.coo1 , multiply(b1.m/(b1.m+b2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
        b1.coo0 = add( b1.coo0 , multiply(b1.m/(b1.m+b2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
        b1.coom1 = add( b1.coom1 , multiply(b1.m/(b1.m+b2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
        b2.coo1 = subtract( b2.coo1 , multiply(b2.m/(b1.m+b2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
        b2.coo0 = subtract( b2.coo0 , multiply(b2.m/(b1.m+b2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
        b2.coom1 = add( b2.coom1 , multiply(b2.m/(b1.m+b2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))


        # if b1 is obj[15] and b2 is obj[28]:
        #     print("v1 ball1:   "+str(vb1))
        #     print("v1 ball2:   "+str(vb2))
        # if b1 is obj[28] and b2 is obj[29]:
        #     #print("impulse")
        #     #print(vb1)
        #     print("cosangle")
        #     print(unitvec_n[0])
        #     print("sinangle")
        #     print(unitvec_n[1])

def Collision_Ball_FixedLine(ball,line):
    # Coordinates of intersection point
    coocol = ClosestPointOnLineSegment(ball.coo1[0],ball.coo1[1],line)
    # Distance between object and intersection point
    distint = sqrt( (coocol[0]-ball.coo1[0])**2 + (coocol[1]-ball.coo1[1])**2 ) # distance between point and fixline

    if distint < ball.r+line.hhalf:
        unitvec_n = divide( array([ ball.coo1[0]-coocol[0] , ball.coo1[1]-coocol[1] , 0 ]) , distint ) # Unit normal vector
        unitvec_t = Rotzvec3_90degcw(unitvec_n) # Unit tangent vector

        coocol_hcomp = add( coocol, multiply( line.hhalf, [ unitvec_n[0], unitvec_n[1] ] ) )

        Contact_dyn_statline(ball,line,coocol_hcomp,unitvec_n,unitvec_t)

        # Old. DO NOT REMOVE. 
        #ball.F_other = add( ball.F_other , multiply( nvec , (2*ball.m/(dt**2))*(ball.r-distint) ) )

        ball.coo1 = ball.coo1 + multiply([unitvec_n[0],unitvec_n[1]],(ball.r+line.hhalf-distint))
        ball.coo0 = ball.coo0 + multiply([unitvec_n[0],unitvec_n[1]],(ball.r+line.hhalf-distint))
        ball.coom1 = ball.coom1 + multiply([unitvec_n[0],unitvec_n[1]],(ball.r+line.hhalf-distint))

def Collision_Ball_Line(ball,line):
    # Coordinates of intersection point
    coocol = ClosestPointOnLineSegment(ball.coo1[0],ball.coo1[1],line)
    # Distance between object and intersection point
    distint = sqrt( (coocol[0]-ball.coo1[0])**2 + (coocol[1]-ball.coo1[1])**2 ) # distance between point and fixline

    if distint < ball.r+line.hhalf:
        unitvec_n = divide( array([ ball.coo1[0]-coocol[0] , ball.coo1[1]-coocol[1] , 0 ]) , distint ) # Unit normal vector
        unitvec_t = Rotzvec3_90degcw(unitvec_n) # Unit tangent vector
    
        coocol_hcomp = add( coocol, multiply( line.hhalf, [ unitvec_n[0], unitvec_n[1] ] ) )

        Contact_dyn_line(ball,line,coocol_hcomp,unitvec_n,unitvec_t) #FELFELFELFELFELFEL FIXA!!!

        pendist = -(distint-(ball.r+line.hhalf))

        ball.coo1 = add( ball.coo1 , multiply(ball.m/(ball.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
        ball.coo0 = add( ball.coo0 , multiply(ball.m/(ball.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
        ball.coom1 = add( ball.coom1 , multiply(ball.m/(ball.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
        line.coo1 = subtract( line.coo1 , multiply(line.m/(ball.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
        line.coo0 = subtract( line.coo0 , multiply(line.m/(ball.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
        line.coom1 = subtract( line.coom1 , multiply(line.m/(ball.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
        line.cooA = subtract( line.cooA , multiply(line.m/(ball.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
        line.cooB = subtract( line.cooB , multiply(line.m/(ball.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))

def Collision_Line_Line(l1,l2):
    coocol_l1A = ClosestPointOnLineSegmentEdgeIndicator(l1.cooA[0],l1.cooA[1],l2)
    coocol_l1B = ClosestPointOnLineSegmentEdgeIndicator(l1.cooB[0],l1.cooB[1],l2)
    coocol_l2A = ClosestPointOnLineSegmentEdgeIndicator(l2.cooA[0],l2.cooA[1],l1)
    coocol_l2B = ClosestPointOnLineSegmentEdgeIndicator(l2.cooB[0],l2.cooB[1],l1)
    distint_l1A = sqrt( (coocol_l1A[0]-l1.cooA[0])**2 + (coocol_l1A[1]-l1.cooA[1])**2 )
    distint_l1B = sqrt( (coocol_l1B[0]-l1.cooB[0])**2 + (coocol_l1B[1]-l1.cooB[1])**2 )
    distint_l2A = sqrt( (coocol_l2A[0]-l2.cooA[0])**2 + (coocol_l2A[1]-l2.cooA[1])**2 )
    distint_l2B = sqrt( (coocol_l2B[0]-l2.cooB[0])**2 + (coocol_l2B[1]-l2.cooB[1])**2 )

    r_pinball = l1.hhalf + l2.hhalf
    fac = 1

    if sqrt( (l1.coo0[0]-l2.coo0[0])**2 + (l1.coo0[1]-l2.coo0[1])**2 ) <= (l1.r + l2.r):

        if distint_l1A <= r_pinball and coocol_l1A[2]==0:
            unitvec_n = divide( [ l1.cooA[0] - coocol_l1A[0] , l1.cooA[1] - coocol_l1A[1] , 0 ] , distint_l1A ) # pekar alltid mot linje 1
            unitvec_t = Rotzvec3_90degcw(unitvec_n)

            coocol_hcomp = add( coocol_l1A, multiply( l2.hhalf, unitvec_n ) )

            Contact_line_line(l1,l2,coocol_hcomp,unitvec_n,unitvec_t)
    
            pendist = (r_pinball-distint_l1A)
            l1.coo1 = add( l1.coo1 , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l1.coo0 = add( l1.coo0 , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l1.coom1 = add( l1.coom1 , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l1.cooA = add( l1.cooA , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l1.cooB = add( l1.cooB , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
    
        if distint_l1B <= r_pinball and coocol_l1B[2]==0:
            unitvec_n = divide( [ l1.cooB[0] - coocol_l1B[0] , l1.cooB[1] - coocol_l1B[1] , 0 ] , distint_l1B ) # pekar alltid mot linje 1
            unitvec_t = Rotzvec3_90degcw(unitvec_n)

            coocol_hcomp = add( coocol_l1B, multiply( l2.hhalf, unitvec_n ) )

            Contact_line_line(l1,l2,coocol_hcomp,unitvec_n,unitvec_t)
        
            pendist = (r_pinball-distint_l1B)
            l1.coo1 = add( l1.coo1 , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l1.coo0 = add( l1.coo0 , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l1.coom1 = add( l1.coom1 , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l1.cooA = add( l1.cooA , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l1.cooB = add( l1.cooB , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
        
        if distint_l2A <= r_pinball and coocol_l2A[2]==0:
            unitvec_n = divide( [ coocol_l2A[0] - l2.cooA[0] , coocol_l2A[1] - l2.cooA[1] , 0 ] , distint_l2A ) # pekar alltid mot linje 1
            unitvec_t = Rotzvec3_90degcw(unitvec_n)

            coocol_hcomp = subtract( coocol_l2A, multiply( l1.hhalf, unitvec_n ) )
            
            Contact_line_line(l1,l2,coocol_hcomp,unitvec_n,unitvec_t)

            pendist = (r_pinball-distint_l2A)
            l1.coo1 = add( l1.coo1 , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l1.coo0 = add( l1.coo0 , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l1.coom1 = add( l1.coom1 , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l1.cooA = add( l1.cooA , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l1.cooB = add( l1.cooB , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
    
        if distint_l2B <= r_pinball and coocol_l2B[2]==0:
            unitvec_n = divide( [ coocol_l2B[0] - l2.cooB[0] , coocol_l2B[1] - l2.cooB[1] , 0 ] , distint_l2B ) # pekar alltid mot linje 1
            unitvec_t = Rotzvec3_90degcw(unitvec_n)

            coocol_hcomp = subtract( coocol_l2B, multiply( l1.hhalf, unitvec_n ) )

            Contact_line_line(l1,l2,coocol_hcomp,unitvec_n,unitvec_t)
            
            pendist = (r_pinball-distint_l2B)
            l1.coo1 = add( l1.coo1 , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l1.coo0 = add( l1.coo0 , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l1.coom1 = add( l1.coom1 , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l1.cooA = add( l1.cooA , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l1.cooB = add( l1.cooB , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))

def Collision_Line_FixedLine(l1,l2):
    coocol_l1A = ClosestPointOnLineSegmentEdgeIndicator(l1.cooA[0],l1.cooA[1],l2)
    coocol_l1B = ClosestPointOnLineSegmentEdgeIndicator(l1.cooB[0],l1.cooB[1],l2)
    coocol_l2A = ClosestPointOnLineSegmentEdgeIndicator(l2.cooA[0],l2.cooA[1],l1)
    coocol_l2B = ClosestPointOnLineSegmentEdgeIndicator(l2.cooB[0],l2.cooB[1],l1)
    distint_l1A = sqrt( (coocol_l1A[0]-l1.cooA[0])**2 + (coocol_l1A[1]-l1.cooA[1])**2 )
    distint_l1B = sqrt( (coocol_l1B[0]-l1.cooB[0])**2 + (coocol_l1B[1]-l1.cooB[1])**2 )
    distint_l2A = sqrt( (coocol_l2A[0]-l2.cooA[0])**2 + (coocol_l2A[1]-l2.cooA[1])**2 )
    distint_l2B = sqrt( (coocol_l2B[0]-l2.cooB[0])**2 + (coocol_l2B[1]-l2.cooB[1])**2 )

    r_pinball = l1.hhalf + l2.hhalf
    fac = 1

    #global colcounter

    if sqrt( (l1.coo0[0]-l2.coo0[0])**2 + (l1.coo0[1]-l2.coo0[1])**2 ) <= (l1.r + l2.r):

        if distint_l1A <= r_pinball and coocol_l1A[2]==0:
            unitvec_n = divide( [ l1.cooA[0] - coocol_l1A[0] , l1.cooA[1] - coocol_l1A[1] , 0 ] , distint_l1A )
            unitvec_t = Rotzvec3_90degcw(unitvec_n)

            coocol_hcomp = add( coocol_l1A, multiply( l2.hhalf, unitvec_n ) )

            Contact_dyn_statline(l1,l2,coocol_hcomp,unitvec_n,unitvec_t)
    
            l1.coo1 = l1.coo1 + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_l1A))
            l1.coo0 = l1.coo0 + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_l1A))
            l1.coom1 = l1.coom1 + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_l1A))
            l1.cooA = l1.cooA + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_l1A))
            l1.cooB = l1.cooB + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_l1A))

            #colcounter = colcounter + 1
            #print(str(colcounter)+" - l1A - "+str(coocol_l1A[2]))
    
        if distint_l1B <= r_pinball and coocol_l1B[2]==0:
            unitvec_n = divide( [ l1.cooB[0] - coocol_l1B[0] , l1.cooB[1] - coocol_l1B[1] , 0 ] , distint_l1B )
            unitvec_t = Rotzvec3_90degcw(unitvec_n)

            coocol_hcomp = add( coocol_l1B, multiply( l2.hhalf, unitvec_n ) )

            Contact_dyn_statline(l1,l2,coocol_hcomp,unitvec_n,unitvec_t)
        
            l1.coo1 = l1.coo1 + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_l1B))
            l1.coo0 = l1.coo0 + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_l1B))
            l1.coom1 = l1.coom1 + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_l1B))
            l1.cooA = l1.cooA + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_l1B))
            l1.cooB = l1.cooB + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_l1B))

            #colcounter = colcounter + 1
            #print(str(colcounter)+" - l1B - "+str(coocol_l1B[2]))
        
        if distint_l2A <= r_pinball and coocol_l2A[2]==0:
            unitvec_n = divide( [ coocol_l2A[0] - l2.cooA[0] , coocol_l2A[1] - l2.cooA[1] , 0 ] , distint_l2A )
            unitvec_t = Rotzvec3_90degcw(unitvec_n)

            coocol_hcomp = subtract( coocol_l2A, multiply( l1.hhalf, unitvec_n ) )
            
            Contact_dyn_statline(l1,l2,coocol_hcomp,unitvec_n,unitvec_t)

            l1.coo1 = l1.coo1 + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_l2A))
            l1.coo0 = l1.coo0 + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_l2A))
            l1.coom1 = l1.coom1 + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_l2A))
            l1.cooA = l1.cooA + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_l2A))
            l1.cooB = l1.cooB + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_l2A))

            #colcounter = colcounter + 1
            #print(str(colcounter)+" - l2A - "+str(coocol_l2A[2]))
    
        if distint_l2B <= r_pinball and coocol_l2B[2]==0:
            unitvec_n = divide( [ coocol_l2B[0] - l2.cooB[0] , coocol_l2B[1] - l2.cooB[1] , 0 ] , distint_l2B )
            unitvec_t = Rotzvec3_90degcw(unitvec_n)

            coocol_hcomp = subtract( coocol_l2B, multiply( l1.hhalf, unitvec_n ) )

            Contact_dyn_statline(l1,l2,coocol_hcomp,unitvec_n,unitvec_t)
            
            l1.coo1 = l1.coo1 + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_l2B))
            l1.coo0 = l1.coo0 + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_l2B))
            l1.coom1 = l1.coom1 + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_l2B))
            l1.cooA = l1.cooA + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_l2B))
            l1.cooB = l1.cooB + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_l2B))

            #colcounter = colcounter + 1
            #print(str(colcounter)+" - l2B - "+str(coocol_l2B[2]))
            
        if obj[10] is l2:
            #print(coocol_l1A)
            Dp.moveto( coocol_l1A[0],coocol_l1A[1] )

class Constraint_SpringDamper:
    def __init__(self,A,B,len0,k,damping):
        self.A = A
        self.B = B
        self.len0 = len0

        self.AB = sqrt(   (self.A.coo1[0]-self.B.coo1[0])**2   +   (self.A.coo1[1]-self.B.coo1[1])**2   )
        self.cosangle = (   (self.B.coo1[0]-self.A.coo1[0])/self.AB   )
        self.sinangle = (   (self.B.coo1[1]-self.A.coo1[1])/self.AB   )

        self.eps = (self.AB-self.len0)/self.len0

        self.vrel = sqrt( (self.A.v0[0] - self.B.v0[0])**2 + (self.A.v0[1] - self.B.v0[1])**2 )
        self.coodist = [abs(self.B.coo1[0] - self.A.coo1[0]), abs(self.B.coo1[1] - self.A.coo1[1])]
        self.line = canvas_1.create_line(self.A.coo1[0],self.A.coo1[1],self.B.coo1[0],self.B.coo1[1])

        self.damping = damping
        self.F_damping = [self.damping*self.cosangle*self.vrel, self.damping*self.sinangle*self.vrel]        
        
        self.len = [self.cosangle*(self.AB-self.len0), self.sinangle*(self.AB-self.len0)]
        self.k = k
        self.update()

    def update(self):
        self.ABm1 = self.AB
        self.AB = sqrt(   (self.A.coo1[0]-self.B.coo1[0])**2   +   (self.A.coo1[1]-self.B.coo1[1])**2   )
        self.deltaAB = self.AB - self.ABm1
        self.cosangle = (   (self.B.coo1[0]-self.A.coo1[0])/self.AB   )
        self.sinangle = (   (self.B.coo1[1]-self.A.coo1[1])/self.AB   )

        self.eps = (self.AB-self.len0)/self.len0

        self.vrel = sqrt( (self.A.v0[0] - self.B.v0[0])**2 + (self.A.v0[1] - self.B.v0[1])**2 )
        self.coodist = array([abs(self.B.coo1[0] - self.A.coo1[0]), abs(self.B.coo1[1] - self.A.coo1[1])])
        self.len = array([self.cosangle*(self.AB-self.len0), self.sinangle*(self.AB-self.len0)])

        self.F_damping = array([self.damping*self.cosangle*self.deltaAB/dt, self.damping*self.sinangle*self.deltaAB/dt])

        self.A.F_other = add( self.A.F_other , array([(self.k*self.len[0] + self.F_damping[0]), (self.k*self.len[1] + self.F_damping[1])]) )
        self.B.F_other = add( self.B.F_other , array([-(self.k*self.len[0] + self.F_damping[0]), -(self.k*self.len[1] + self.F_damping[1])]) )

        # Draw
        canvas_1.coords(self.line,self.A.coo1[0],self.A.coo1[1],self.B.coo1[0],self.B.coo1[1])
        canvas_1.itemconfigure(self.line,fill="black",width=2)
        if self.eps > 0.15:
            canvas_1.itemconfigure(self.line,fill="red")
        else:
            if self.eps < -0.15:
                canvas_1.itemconfigure(self.line,fill="blue",width=(1-self.eps)**6)
            
class Constraint_Distance:
    def __init__(self,A,B,len0):
        self.A = A
        self.B = B
        self.len0 = len0

        self.AB = sqrt(   (self.A.coo1[0]-self.B.coo1[0])**2   +   (self.A.coo1[1]-self.B.coo1[1])**2   )
        self.cosangle = (   (self.B.coo1[0]-self.A.coo1[0])/self.AB   )
        self.sinangle = (   (self.B.coo1[1]-self.A.coo1[1])/self.AB   )

        self.eps = (self.AB-self.len0)/self.len0

        self.line = canvas_1.create_line(self.A.coo1[0],self.A.coo1[1],self.B.coo1[0],self.B.coo1[1])

        self.len = [self.cosangle*(self.AB-self.len0), self.sinangle*(self.AB-self.len0)]
        self.update()

    def update(self):
        self.ABm1 = self.AB
        self.AB = sqrt(   (self.A.coo1[0]-self.B.coo1[0])**2   +   (self.A.coo1[1]-self.B.coo1[1])**2   )
        self.cosangle = (   (self.B.coo1[0]-self.A.coo1[0])/self.AB   )
        self.sinangle = (   (self.B.coo1[1]-self.A.coo1[1])/self.AB   )

        self.eps = (self.AB-self.len0)/self.len0

        self.len = array([self.cosangle*(self.AB-self.len0), self.sinangle*(self.AB-self.len0)])

        self.A.coo1 = [ self.A.coo1[0] + self.len[0]/2 , self.A.coo1[1] + self.len[1]/2 ]
        #self.A.coo0 = [ self.A.coo0[0] + self.len[0]/2 , self.A.coo0[1] + self.len[1]/2 ]
        #self.A.coom1 = [ self.A.coom1[0] + self.len[0]/2 , self.A.coom1[1] + self.len[1]/2 ]
        self.B.coo1 = [ self.B.coo1[0] - self.len[0]/2 , self.B.coo1[1] - self.len[1]/2 ]
        #self.B.coo0 = [ self.B.coo0[0] - self.len[0]/2 , self.B.coo0[1] - self.len[1]/2 ]
        #self.B.coom1 = [ self.B.coom1[0] - self.len[0]/2 , self.B.coom1[1] - self.len[1]/2 ]
    

        # Draw
        canvas_1.coords(self.line,self.A.coo1[0],self.A.coo1[1],self.B.coo1[0],self.B.coo1[1])
        canvas_1.itemconfigure(self.line,fill="black",width=2)
        if self.eps > 0.15:
            canvas_1.itemconfigure(self.line,fill="red")
        else:
            if self.eps < -0.15:
                canvas_1.itemconfigure(self.line,fill="blue",width=(1-self.eps)**6)

class DebugPoint:
    def __init__(self,x,y):
        self.coo1 = array([x, y])
        self.r = 8
        self.ball = canvas_1.create_oval(x-self.r,y-self.r,x+self.r,y+self.r,fill="red",outline="white")
    def moveto(self,x,y):
        self.coo1 = array([x, y])
        canvas_1.coords(self.ball,self.coo1[0]-self.r,self.coo1[1]-self.r,self.coo1[0]+self.r,self.coo1[1]+self.r)

class Object_FixedPoint:
    def __init__(self,canvas,x,y):
        self.coofix = [x,y]
        self.coo1 = array([x, y])
        self.coo0 = array([x, y])
        self.coom1 = array([x, y])
        self.r = 4
        self.v0 = array([0,0])
        self.F_other = array([0,0])
        self.ball = canvas.create_oval(x-self.r,y-self.r,x+self.r,y+self.r,fill="green")
    def integrate(self):
        self.coo1 = self.coofix
        return

class Object_Ball:
    def __init__(self,canvas,x,y,d,rho,C_el,C_fric):
        # Indata
        self.coo1 = array([x, y])
        self.coo0 = self.coo1
        self.coom1 = self.coo1
        self.theta1 = 0.0
        self.theta0 = 0.0
        self.thetam1 = 0.0
        self.d = d
        self.r = d/2
        self.A = 0.25*pi*d**2
        self.V = 0.16666*pi*d**3
        self.rho = rho
        self.C_fric = C_fric
        self.C_el = C_el
        self.v0 = [0.0, 0.0]
        self.C_D = 0.47
        self.coothetaline = add( self.coo1, Rotvec2([self.r,0],self.theta1) )
        #self.dampen = 0

        self.vel = subtract(self.coo1,self.coo0)/dt
        self.absvel = sqrt( (self.vel[0])**2 + (self.vel[1])**2 )
        self.unitvec_vel = [Safediv(self.vel[0],self.absvel),Safediv(self.vel[1],self.absvel)]
        self.A_ABfacingdirection = pi*self.r**2

        # Inertia
        self.m = self.rho*self.V
        self.invm = 1/self.m
        self.I = 0.16667*self.m*self.d*self.d
        self.invImat = [ [1/self.I,0,0] , [0,1/self.I,0] , [0,0,1/self.I] ]
        # Forces
        self.F_other = array([0.0,0.0])
        self.F_other_abs = round(sqrt(self.F_other[0]**2 + self.F_other[1]**2)/10000)
        self.tau_other = 0.0
        self.F_g = array([0.0, 0.0])
        self.F_D = array([0.0, 0.0])
        self.F_R = array([0.0, 0.0])
        # Draw
        self.canvas = canvas
        self.ball = canvas.create_oval(x-self.r,y-self.r,x+self.r,y+self.r,fill="white")
        self.thetaline = canvas.create_line(self.coo1[0],self.coo1[1],self.coothetaline[0],self.coothetaline[1])
        #self.arrow_F = canvas.create_line(self.coo1[0],self.coo1[1],self.F_other[0],self.F_other[1],arrow=LAST,fill="grey")
        #self.text_F = canvas_1.create_text(self.coo0[0],self.coo0[1],text=str(self.F_other_abs),font=("arial",8),fill="grey")
        
        self.arrow_g = [ sf_farrow1*tanh((self.m*grav[0])/sf_farrow2) , sf_farrow1*tanh((self.m*grav[1])/sf_farrow2) ]
        self.arrow_F_D = [ (sf_farrow1*tanh(self.F_D[0]/sf_farrow2)) , (sf_farrow1*tanh(self.F_D[1]/sf_farrow2)) ]
        self.arrow_g_draw = canvas.create_line(self.coo1[0],self.coo1[1],self.coo1[0]+self.arrow_g[0],self.coo1[1]+self.arrow_g[1],arrow=LAST,fill="red")
        #self.offset_text_coo = -20
        #self.text_coo = canvas_1.create_text(self.coo0[0],self.coo0[1]+self.offset_text_coo,text="("+str(round(self.coo0[0]))+", "+str(round(self.coo0[1]))+")",font=("arial",8))        
        self.arrow_F_D_draw = canvas.create_line(self.coo1[0],self.coo1[1],self.coo1[0]-self.arrow_F_D[0],self.coo1[1]-self.arrow_F_D[1],arrow=LAST,fill="blue")
        

    def colcheck(self):
        if self.coo1[1] + self.r >= h_cnv:
            self.F_R = -100000000*((self.coo1[1] + self.r) - h_cnv)
            self.coo1[1] = self.coo1[1] + 0.5*(self.F_R/self.m)*dt**2 #0.5 här eller ej?
            canvas_1.itemconfigure(self.ball,fill="black")
        else:
            if self.coo1[0] + self.r >= w_cnv:
                self.F_R = -100000000*((self.coo1[0] + self.r) - w_cnv)
                self.coo1[0] = self.coo1[0] + 0.5*(self.F_R/self.m)*dt**2 #0.5 här eller ej?
                canvas_1.itemconfigure(self.ball,fill="black")
            else:
                if self.coo1[0] - self.r <= 0:
                    self.F_R = -100000000*((self.coo1[0] - self.r) - 0)
                    self.coo1[0] = self.coo1[0] + 0.5*(self.F_R/self.m)*dt**2 #0.5 här eller ej?
                    canvas_1.itemconfigure(self.ball,fill="black")
                else:
                    canvas_1.itemconfigure(self.ball,fill="white")

    def draw(self):
        
        self.canvas.coords(self.ball,self.coo1[0]-self.r,self.coo1[1]-self.r,self.coo1[0]+self.r,self.coo1[1]+self.r)
        self.canvas.coords(self.thetaline,self.coo1[0],self.coo1[1],self.coothetaline[0],self.coothetaline[1])

        #self.canvas.coords(self.arrow_F,self.coo1[0],self.coo1[1],self.coo1[0]+Safediv(self.F_other[0],self.F_other_abs)*viz_sf,self.coo1[1]+Safediv(self.F_other[1],self.F_other_abs)*viz_sf)
        #self.canvas.coords(self.text_F,self.coo1[0]+Safediv(self.F_other[0],self.F_other_abs)*viz_sf,self.coo1[1]+Safediv(self.F_other[1],self.F_other_abs)*viz_sf)

        #self.canvas.itemconfigure(self.text_F,text=str(round(self.F_other_abs/10000)))
        
        
        self.arrow_g = [ sf_farrow1*tanh((self.m*grav[0])/sf_farrow2) , sf_farrow1*tanh((self.m*grav[1])/sf_farrow2) ]
        self.arrow_F_D = [ (sf_farrow1*tanh(self.F_D[0]/sf_farrow2)) , (sf_farrow1*tanh(self.F_D[1]/sf_farrow2)) ]
        self.canvas.coords(self.arrow_g_draw,self.coo1[0],self.coo1[1],self.coo1[0]+self.arrow_g[0],self.coo1[1]+self.arrow_g[1])
        self.canvas.coords(self.arrow_F_D_draw,self.coo1[0],self.coo1[1],self.coo1[0]-self.arrow_F_D[0],self.coo1[1]-self.arrow_F_D[1])
        

        #self.canvas.coords(self.text_coo,self.coo0[0],self.coo0[1]+self.offset_text_coo)
        #self.canvas.itemconfigure(self.text_coo,text="("+str(round(self.coo0[0]))+", "+str(round(self.coo0[1]))+")")

    def integrate(self):

        self.vel = subtract(self.coo1,self.coo0)/dt
        self.absvel = sqrt( (self.vel[0])**2 + (self.vel[1])**2 )
        self.unitvec_vel = [Safediv(self.vel[0],self.absvel),Safediv(self.vel[1],self.absvel)]
        self.A_ABfacingdirection = pi*self.r**2

        # Drag force in the chosen medium (air, water, etc.)
        self.F_D = multiply( 0.5*rho_medium*0.50*self.A_ABfacingdirection , multiply( [ (self.vel[0])**2 , (self.vel[1])**2 ] , self.unitvec_vel ) )
        F_crit = self.m*self.absvel/(2*dt)
        if sqrt( (self.F_D[0])**2 + (self.F_D[1])**2 ) >= F_crit:
            self.F_D = multiply(F_crit,self.unitvec_vel)
            #print("                                      Drag force limit exceeded" + str(self.F_D))
        self.F_other = subtract( self.F_other , self.F_D )
        # LÄGG TILL LUFTMOTSTÅND MOT ROTATION

        #self.coom1 = subtract( self.coo0, multiply(self.v0,dt) )
        self.coom1 = self.coo0
        #print(self.coom1)
        self.coo0 = self.coo1
        #print(self.coo0)
        self.thetam1 = self.theta0
        self.theta0 = self.theta1
        #self.v0 = [0,0]#divide(subtract(self.coo0,self.coo1), dt)
        self.F_g = multiply(self.m,grav)
        self.F_other = add( self.F_other, self.F_g )
        
        self.coo1 = [ 2*self.coo0[0] - self.coom1[0] + (self.F_other[0]/self.m)*dt**2 , 2*self.coo0[1] - self.coom1[1] + (self.F_other[1]/self.m)*dt**2 ]
        self.theta1 = 2*self.theta0 - self.thetam1 + (self.tau_other/self.I)*dt**2
        self.coothetaline = add( self.coo1, Rotvec2([self.r,0],self.theta1) )


        self.F_other_abs = sqrt(self.F_other[0]**2 + self.F_other[1]**2)
        #self.coo1 = subtract(self.coo1 , self.dampen*subtract(self.coo1 , self.coo0))

        #self.dampen = 0

        self.colcheck()
        self.draw() # ändra så att draw är en allmän funktion som kallas efter den globala kollisionskollen

        #ändra så att den tar hänsyn till studsen
        #self.v0 = [(self.coo1[0] - self.coo0[0])/dt, (self.coo1[1] - self.coo0[1])/dt]

class Object_Box:
    def __init__(self,canvas,x,y,w,h,t,rho):
        # General
        self.canvas = canvas
        # Geometrical properties
        self.w = w
        self.h = h
        self.t = t
        # Coordinates
        self.coo1 = array([ x, y, x-0.5*w , y-0.5*h , x+0.5*w , y-0.5*h , x+0.5*w , y+0.5*h , x-0.5*w , y+0.5*h ])
        self.k = array([ Safediv((self.coo1[5]-self.coo1[3]),(self.coo1[4]-self.coo1[2])) , Safediv((self.coo1[7]-self.coo1[5]),(self.coo1[6]-self.coo1[4])) , Safediv((self.coo1[9]-self.coo1[7]),(self.coo1[8]-self.coo1[6])) , Safediv((self.coo1[3]-self.coo1[9]),(self.coo1[2]-self.coo1[8])) ])
        #print(self.k)
        #self.coo1 = [x,y]
        self.coo0 = self.coo1
        self.coom1 = self.coo1
        self.theta1 = 0
        self.theta0 = self.theta1
        self.thetam1 = self.theta1
        # Forces and torques
        self.F_g = [0.0,0.0]
        self.F_other = [0.0,0.0]
        self.tau_other = 0.0
        # Inertia
        self.m = w*h*t*rho
        self.I = 0.08333*self.m*(w*w+h*h)
        # Draw
        self.lineAB = canvas.create_line(self.coo0[2],self.coo0[3],self.coo0[4],self.coo0[5])
        self.lineBC = canvas.create_line(self.coo0[4],self.coo0[5],self.coo0[6],self.coo0[7])
        self.lineCD = canvas.create_line(self.coo0[6],self.coo0[7],self.coo0[8],self.coo0[9])
        self.lineDA = canvas.create_line(self.coo0[8],self.coo0[9],self.coo0[2],self.coo0[3])
        self.ball = canvas.create_oval(x-5,y-5,x+5,y+5,fill="white")


    def draw(self):
        self.canvas.coords( self.lineAB,self.coo1[2],self.coo1[3],self.coo1[4],self.coo1[5] )
        self.canvas.coords( self.lineBC,self.coo1[4],self.coo1[5],self.coo1[6],self.coo1[7] )
        self.canvas.coords( self.lineCD,self.coo1[6],self.coo1[7],self.coo1[8],self.coo1[9] )
        self.canvas.coords( self.lineDA,self.coo1[8],self.coo1[9],self.coo1[2],self.coo1[3] )
        #self.canvas.coords( self.ball,self.coo1[0]-5,self.coo1[1]-5,self.coo1[0]+5,self.coo1[1]+5 )
        self.canvas.coords( self.ball,self.coo1[0]-5,self.coo1[1]-5,self.coo1[0]+5,self.coo1[1]+5 )


    def integrate(self):

        self.F_g = multiply(self.m,grav)
        self.F_other = add( self.F_other, self.F_g )

        self.coom1 = self.coo0
        self.coo0 = self.coo1
        self.v0 = array([0,0])#divide(subtract(self.coo0,self.coo1), dt)
        self.coo1 = [ 2*self.coo0[0] - self.coom1[0] + (self.F_other[0]/self.m)*dt**2 , 2*self.coo0[1] - self.coom1[1] + (self.F_other[1]/self.m)*dt**2 ]
        self.coo1 = [ self.coo1[0], self.coo1[1], self.coo1[0]-0.5*self.w, self.coo1[1]-0.5*self.h, self.coo1[0]+0.5*self.w, self.coo1[1]-0.5*self.h, self.coo1[0]+0.5*self.w, self.coo1[1]+0.5*self.h, self.coo1[0]-0.5*self.w, self.coo1[1]+0.5*self.h ]
        #self.coo1 = [ 0,0,0,0,0,0,0,0, 2*self.coo0[8] - self.coom1[8] + (self.F_other[0]/self.m)*dt**2 , 2*self.coo0[9] - self.coom1[9] + (self.F_other[1]/self.m + grav)*dt**2 ]
        #f1 = (2*self.coo0[0] - self.coom1[0] + (self.F_other[0]/self.m)*dt**2)
        #f2 = (2*self.coo0[1] - self.coom1[1] + (self.F_other[1]/self.m + grav)*dt**2)
        #self.coo1[0] = f1
        #self.coo1[1] = f2
        
        #self.coom1 = self.coo0
        #print(self.coom1)
        #self.coo0 = self.coo1
        #print(self.coo0)
        #self.thetam1 = self.theta0
        #self.theta0 = self.theta1
        #self.coo1[8] = (2*self.coo0[8] - self.coom1[8] + (self.F_other[0]/self.m)*dt**2)
        #self.coo1[9] = (2*self.coo0[9] - self.coom1[9] + (self.F_other[1]/self.m + grav)*dt**2)
        #self.coo1[9] = 2*self.coo0[9] - self.coom1[9] + (self.F_other[1]/self.m + grav)*dt**2
        # self.coo1[0] = self.coo1[8]-0.5*self.w
        # self.coo1[1] = self.coo1[9]-0.5*self.h
        # self.coo1[2] = self.coo1[8]+0.5*self.w
        # self.coo1[3] = self.coo1[9]-0.5*self.h
        # self.coo1[4] = self.coo1[8]+0.5*self.w
        # self.coo1[5] = self.coo1[9]+0.5*self.h
        # self.coo1[6] = self.coo1[8]-0.5*self.w
        # self.coo1[7] = self.coo1[9]+0.5*self.h
        
       

        #self.theta1 = 2*self.theta0 - self.thetam1 + (self.tau_other/self.I)*dt**2
        #print(len(self.coo0))
        #print(self.coo1)
        #print(self.coo0[9]-self.coom1[9])
        #print(self.F_other)

        self.draw()

class Object_FixedLine:
    def __init__(self,canvas,xA,yA,xB,yB,h,C_el,C_fric):
        # Indata
        self.C_fric = C_fric
        self.C_el = C_el

        self.h = h
        self.hhalf = h/2

        self.F_other = [0,0]
        self.F_R = array([0.0, 0.0])
        self.cooA = [xA, yA]
        self.cooB = [xB, yB]
        self.coo0 = [(xA+xB)/2, (yA+yB)/2]
        # Bättre att beräkna dessa en gång, istället för varje gång under kontaktsökningen
        self.coomax = [ max(self.cooA[0],self.cooB[0]) , max(self.cooA[1],self.cooB[1]) ]
        self.coomin = [ min(self.cooA[0],self.cooB[0]) , min(self.cooA[1],self.cooB[1]) ]

        self.AB = sqrt( (xB-xA)**2 + (yB-yA)**2 )
        self.r = self.AB/2
        self.ang = arccos( (xB-xA)/self.AB )
        #self.coomidp = [ self.cooA[0] + 0.5*(self.cooB[0]-self.cooA[0]) , self.cooA[1] + 0.5*(self.cooB[1]-self.cooA[1]) ]
        self.unitvec_n = Rotvec2( [1,0] , self.ang-pi/2 )
        #self.normalvec = add( Rotvec2( [20,0] , self.ang-pi/2 ) , self.coomidp )
        #self.nvec = [0,0]
        self.k = Safediv( (self.cooB[1]-self.cooA[1]) , (self.cooB[0]-self.cooA[0]) )
        self.coocol = [0,0]
        self.distint = 0
        #self.text_distint = canvas_1.create_text(0,0,text="",font=("arial",8))

        # Forces
        self.F_other = [0.0,0.0]
        # Draw
        self.canvas = canvas
        self.lineAB = canvas.create_line(self.cooA[0],self.cooA[1],self.cooB[0],self.cooB[1],width=h)
        self.ballA = canvas.create_oval(self.cooA[0]-self.hhalf,self.cooA[1]-self.hhalf,self.cooA[0]+self.hhalf,self.cooA[1]+self.hhalf,outline="black",fill="black")
        self.ballB = canvas.create_oval(self.cooB[0]-self.hhalf,self.cooB[1]-self.hhalf,self.cooB[0]+self.hhalf,self.cooB[1]+self.hhalf,outline="black",fill="black")
        #self.lineN = canvas.create_line(self.coomidp[0],self.coomidp[1],self.normalvec[0],self.normalvec[1],arrow=LAST)
        #self.linecol = canvas.create_line( 0,0,1,0, dash = (2,2) )


    def integrate(self):
        pass

    # def Distto(self,coop): # Intersection between normal line and main line, where normal line is at its closest point to object

    #     # Distance
    #     self.coocol[0] = ( coop[1]+coop[0]/self.k-self.cooA[1]+self.cooA[0]*self.k )/(self.k+1/self.k) # 2 line intersection
    #     self.coocol[0] = ( (self.coocol[0]>=self.coomin[0]) and (self.coocol[0]<=self.coomax[0]) )*self.coocol[0] + ( (self.coocol[0]<=self.coomin[0]) )*self.coomin[0] + ( (self.coocol[0]>=self.coomax[0]) )*self.coomax[0]
    #     self.coocol[1] = self.k*(self.coocol[0]-self.cooA[0])+self.cooA[1] # intersection x into fixline equation
    #     self.distint = sqrt( (self.coocol[0]-coop[0])**2 + (self.coocol[1]-coop[1])**2 ) # distance between point and fixline

    #     # Normal vector
    #     self.nvec = divide( [ coop[0]-self.coocol[0] , coop[1]-self.coocol[1] ] , sqrt( (coop[0]-self.coocol[0])**2 + (coop[1]-self.coocol[1])**2 ) )
        
    #     #print("coocol: "+str(self.coocol)+ "          objcoo: "+str(coop)+ "          dist: "+str(self.distint))

    #     # Draw
    #     self.canvas.coords(self.linecol,self.coocol[0],self.coocol[1],coop[0],coop[1])
    #     self.canvas.coords(self.text_distint, self.coocol[0] + 0.5*(coop[0] - self.coocol[0]) , self.coocol[1] + 0.5*(coop[1] - self.coocol[1]) )
    #     self.canvas.itemconfigure(self.text_distint,text=round(self.distint))

    #     #print(self.coocol)
        
    #     # Return value
    #     return self.distint    

class Object_Line:
    def __init__(self,canvas,xA,yA,xB,yB,h,t,rho,C_el,C_fric):
        # Indata
        self.canvas = canvas
        self.AB = sqrt( (xB-xA)**2 + (yB-yA)**2 )
        self.vec_AB = [ (xB-xA) , (yB-yA) ]
        self.r = self.AB/2
        self.C_fric = C_fric
        self.C_el = C_el

        self.h = h
        self.t = t
        self.hhalf = h/2

        self.coo1 = array([xA + 0.5*(xB-xA), yA + 0.5*(yB-yA)])
        self.coo0 = self.coo1
        self.coom1 = self.coo1
        self.theta1 = arcsin( Safediv( abs(yB-yA) , self.AB ) )
        if xB<xA:
            self.theta1 = pi - self.theta1
        if yB<yA:
            self.theta1 = - self.theta1

        self.vel = subtract(self.coo1,self.coo0)/dt
        self.absvel = sqrt( (self.vel[0])**2 + (self.vel[1])**2 )
        self.unitvec_vel = [Safediv(self.vel[0],self.absvel),Safediv(self.vel[1],self.absvel)]
        self.unitvec_velperp = Rotvec2_90degcw(self.unitvec_vel)
        self.vec_ABfacingdirection = Vecproj_vec2(self.vec_AB,self.unitvec_velperp)
        self.A_ABfacingdirection = self.t*sqrt( (self.vec_ABfacingdirection[0])**2 + (self.vec_ABfacingdirection[1])**2 )


        self.theta0 = self.theta1
        self.thetam1 = self.theta1
        self.F_other = [0,0]
        self.F_D = [0,0]
        self.F_R = array([0.0, 0.0]) # Behövs verkligen denna?
        self.tau_other = 0
        self.cooA = [xA, yA]
        self.cooB = [xB, yB]
        self.cooAloc = subtract(self.cooA,self.coo1)
        self.cooBloc = subtract(self.cooB,self.coo1)
         # Bättre att beräkna dessa en gång, istället för varje gång under kontaktsökningen
        self.coomax = [ max(self.cooA[0],self.cooB[0]) , max(self.cooA[1],self.cooB[1]) ] #ta bort?
        self.coomin = [ min(self.cooA[0],self.cooB[0]) , min(self.cooA[1],self.cooB[1]) ] #ta bort?
        #print(self.cooAloc)
        #print(self.cooA)
        #print(self.cooB)

        self.V = self.AB*h*t
        self.m = rho*self.V
        self.invm = 1/self.m
        self.I = 0.08333*self.m*(h*h+self.AB**2)
        self.invImat = [ [1/self.I,0,0] , [0,0,0] , [0,0,1/self.I] ]
        
        self.ang = arccos( (xB-xA)/self.AB )
        #print(self.theta1*180/pi)

        #if self.cooA[1]>self.cooB[1]:
        #    self.ang = -self.ang

        #print(self.ang*180/pi)
        #self.coomidp = [ self.cooA[0] + 0.5*(self.cooB[0]-self.cooA[0]) , self.cooA[1] + 0.5*(self.cooB[1]-self.cooA[1]) ]
        self.unitvec_n = Rotvec2( [1,0] , self.ang-pi/2 )
        #self.normalvec = add( Rotvec2( [20,0] , self.ang-pi/2 ) , self.coomidp )
        #self.nvec = [0,0]
        self.k = Safediv( (self.cooB[1]-self.cooA[1]) , (self.cooB[0]-self.cooA[0]) )
        self.coocol = [0,0]
        self.distint = 0
        #self.text_distint = canvas_1.create_text(0,0,text="",font=("arial",8))

        # Forces
        self.F_other = [0.0,0.0]
        # Draw
        self.canvas = canvas
        self.lineAB = canvas.create_line(self.cooA[0],self.cooA[1],self.cooB[0],self.cooB[1],width=h,fill="grey")
        self.ballA = canvas.create_oval(self.cooA[0]-self.hhalf,self.cooA[1]-self.hhalf,self.cooA[0]+self.hhalf,self.cooA[1]+self.hhalf,outline="grey",fill="grey")
        self.ballB = canvas.create_oval(self.cooB[0]-self.hhalf,self.cooB[1]-self.hhalf,self.cooB[0]+self.hhalf,self.cooB[1]+self.hhalf,outline="grey",fill="grey")
        
        self.arrow_g = [ sf_farrow1*tanh((self.m*grav[0])/sf_farrow2) , sf_farrow1*tanh((self.m*grav[1])/sf_farrow2) ]
        self.arrow_F_D = [ (sf_farrow1*tanh(self.F_D[0]/sf_farrow2)) , (sf_farrow1*tanh(self.F_D[1]/sf_farrow2)) ]
        self.arrow_g_draw = canvas.create_line(self.coo1[0],self.coo1[1],self.coo1[0]+self.arrow_g[0],self.coo1[1]+self.arrow_g[1],arrow=LAST,fill="red")
        #self.offset_text_coo = -20
        #self.text_coo = canvas_1.create_text(self.coo0[0],self.coo0[1]+self.offset_text_coo,text="("+str(round(self.coo0[0]))+", "+str(round(self.coo0[1]))+")",font=("arial",8))        
        self.arrow_F_D_draw = canvas.create_line(self.coo1[0],self.coo1[1],self.coo1[0]-self.arrow_F_D[0],self.coo1[1]-self.arrow_F_D[1],arrow=LAST,fill="blue")
        
        #self.lineN = canvas.create_line(self.coomidp[0],self.coomidp[1],self.normalvec[0],self.normalvec[1],arrow=LAST)
        #self.linecol = canvas.create_line( 0,0,1,0, dash = (2,2) )
        



    def integrate(self):
        
        self.vec_AB = [ (self.cooB[0]-self.cooA[0]) , (self.cooB[1]-self.cooA[1]) ]
        # Ersätt med en Safedic för vektorer. Nämnaren är 0 vid programstart.
        self.vel = subtract(self.coo1,self.coo0)/dt
        self.absvel = sqrt( (self.vel[0])**2 + (self.vel[1])**2 )
        self.unitvec_vel = [Safediv(self.vel[0],self.absvel),Safediv(self.vel[1],self.absvel)]
        self.unitvec_velperp = Rotvec2_90degcw(self.unitvec_vel)
        self.vec_ABfacingdirection = Vecproj_vec2(self.vec_AB,self.unitvec_velperp)
        self.A_ABfacingdirection = self.t*sqrt( (self.vec_ABfacingdirection[0])**2 + (self.vec_ABfacingdirection[1])**2 )

        # Drag force in the chosen medium (air, water, etc.)
        self.F_D = multiply( 0.5*rho_medium*1.98*self.A_ABfacingdirection , multiply( [ (self.vel[0])**2 , (self.vel[1])**2 ] , self.unitvec_vel ) )
        F_crit = self.m*self.absvel/(2*dt)
        if sqrt( (self.F_D[0])**2 + (self.F_D[1])**2 ) >= F_crit:
            self.F_D = multiply(F_crit,self.unitvec_vel)
            #print("                                      Drag force limit exceeded" + str(self.F_D))
        self.F_other = subtract( self.F_other , self.F_D )
        # LÄGG TILL LUFTMOTSTÅND MOT ROTATION
        
        self.coom1 = self.coo0
        self.coo0 = self.coo1
        self.thetam1 = self.theta0
        self.theta0 = self.theta1
        #print(self.coo0)
        #self.v0 = [0,0]#divide(subtract(self.coo0,self.coo1), dt)

        self.F_g = multiply(self.m,grav)
        self.F_other = add( self.F_other, self.F_g )

        self.coo1 = [ 2*self.coo0[0] - self.coom1[0] + (self.F_other[0]/self.m)*dt**2 , 2*self.coo0[1] - self.coom1[1] + (self.F_other[1]/self.m)*dt**2 ]
        self.theta1 = 2*self.theta0 - self.thetam1 + (self.tau_other/self.I)*dt**2
        self.k = Safediv( (self.cooB[1]-self.cooA[1]) , (self.cooB[0]-self.cooA[0]) )
        self.cooAloc = Rotvec2([-0.5*self.AB,0],self.theta1)
        self.cooBloc = Rotvec2([0.5*self.AB,0],self.theta1)
        self.cooA = add(self.coo1 , self.cooAloc)
        self.cooB = add(self.coo1 , self.cooBloc)
        self.coomax = [ max(self.cooA[0],self.cooB[0]) , max(self.cooA[1],self.cooB[1]) ] #ta bort?
        self.coomin = [ min(self.cooA[0],self.cooB[0]) , min(self.cooA[1],self.cooB[1]) ] #ta bort?

        #print(self.cooA)
        #print(self.cooB)
        self.canvas.coords(self.lineAB,self.cooA[0],self.cooA[1],self.cooB[0],self.cooB[1])
        self.canvas.coords(self.ballA,self.cooA[0]-self.hhalf,self.cooA[1]-self.hhalf,self.cooA[0]+self.hhalf,self.cooA[1]+self.hhalf)
        self.canvas.coords(self.ballB,self.cooB[0]-self.hhalf,self.cooB[1]-self.hhalf,self.cooB[0]+self.hhalf,self.cooB[1]+self.hhalf)

        self.arrow_g = [ sf_farrow1*tanh((self.m*grav[0])/sf_farrow2) , sf_farrow1*tanh((self.m*grav[1])/sf_farrow2) ]
        self.arrow_F_D = [ (sf_farrow1*tanh(self.F_D[0]/sf_farrow2)) , (sf_farrow1*tanh(self.F_D[1]/sf_farrow2)) ]
        self.canvas.coords(self.arrow_g_draw,self.coo1[0],self.coo1[1],self.coo1[0]+self.arrow_g[0],self.coo1[1]+self.arrow_g[1])
        self.canvas.coords(self.arrow_F_D_draw,self.coo1[0],self.coo1[1],self.coo1[0]-self.arrow_F_D[0],self.coo1[1]-self.arrow_F_D[1])

        #print(self.theta1*180/pi)
        #print(abs(cos(self.theta1)))
    
class Object_Trace:
    def __init__(self,canvas,target,freq,nlinemax):
        self.canvas = canvas
        self.target = target
        self.coom1 = target.coom1
        self.freq = freq
        self.nlinemax = nlinemax
        self.counter = 0
        self.linelist = []

    def integrate(self):
        self.counter = self.counter + 1
        if self.counter >= self.freq:
            self.linelist.append( self.canvas.create_line( self.coom1[0],self.coom1[1],self.target.coo1[0],self.target.coo1[1],fill="grey", dash = (2,2) ) )
            if len(self.linelist)>self.nlinemax-1:
                canvas_1.delete(self.linelist[0])
                self.linelist.pop(0)
            self.coom1 = self.target.coom1
            self.counter = 0

class Object_ShowPhysics:
    def __init__(self,canvas,target):
        self.canvas = canvas
        self.target = target
        #self.F = target.F_other
        self.coo0 = target.coo0
        self.coo1 = target.coo1
        self.vel = subtract(self.coo1,self.coo0)/dt
        self.l_arrow = 40
        self.offset_arrow = 30
        self.absvel = sqrt( (self.vel[0])**2 + (self.vel[1])**2 )
        #self.absF = sqrt( (self.F[0])**2 + (self.F[1])**2 )
        self.velarrow = multiply( self.l_arrow , Safediv(self.vel,self.absvel) )
        #self.Farrow = multiply( self.l_arrow , Safediv(self.F,self.absF) )
        self.offset = 65
        self.velarrow2 = canvas.create_line(self.coo1[0],self.coo1[1]-self.offset,self.coo1[0]+self.velarrow[0],self.coo1[1]+self.offset+self.velarrow[1],arrow=LAST,fill="purple")
        self.text_vel = canvas_1.create_text(0,0,text="",font=("arial",8),fill="purple")
        #self.Farrow2 = canvas.create_line(self.coo1[0],self.coo1[1]-self.offset,self.coo1[0]+self.Farrow[0],self.coo1[1]+self.offset+self.Farrow[1],arrow=LAST,fill="grey")


    def integrate(self):
        self.coo0 = self.target.coo0
        self.coo1 = self.target.coo1
        self.theta0 = self.target.theta0
        self.theta1 = self.target.theta1
        self.vel = subtract(self.coo1,self.coo0)/dt
        #self.F = self.target.F_other
        self.omega = subtract(self.theta1,self.theta0)/dt
        self.absvel = sqrt( (self.vel[0])**2 + (self.vel[1])**2 )
        self.velarrow = multiply( self.l_arrow , Safediv(self.vel,self.absvel) )
        self.Ektrans = 0.5*self.target.m*self.absvel**2
        self.Ekrot = 0.5*self.target.I*self.omega**2
        canvas_1.coords(self.velarrow2,self.coo1[0]-self.velarrow[0]/2,self.coo1[1]-self.offset_arrow-self.velarrow[1]/2,self.coo1[0]+self.velarrow[0]/2,self.coo1[1]-self.offset_arrow+self.velarrow[1]/2)
        #canvas_1.coords(self.Farrow2,self.coo1[0],self.coo1[1]-self.offset_arrow,self.coo1[0]+self.Farrow[0],self.coo1[1]-self.offset_arrow+self.Farrow[1])
        #canvas_1.coords(self.Farrow2,self.coo1[0],self.coo1[1]-self.offset_arrow,self.coo1[0]+self.F[0],self.coo1[1]-self.offset_arrow+self.F[1])
        
        canvas_1.itemconfigure(self.text_vel,text="x, y = "+str(round(self.coo0[0]))+", "+str(round(self.coo0[1]))+"\nvel = "+str(round(self.absvel,3))+"\nomega = "+str(round(self.omega,3))+"\nvel_t = "+str(round(self.omega*self.target.r,3))+"\nE_k = "+str(round(self.Ektrans+self.Ekrot)))
        canvas_1.coords(self.text_vel,self.coo1[0]+self.offset,self.coo1[1]-self.offset)

class Object_ShowDistance_Point_LineSegment:
    def __init__(self,o1,o1_coo,o2):
        self.o1 = o1
        self.o1_coo = o1_coo
        self.o1_point = eval("self.o1." + o1_coo)
        #print(eval("self.o1." + o1_coo))
        #print(self.o1_point)
        self.o2 = o2
        self.coocol = [0,0]
        self.distint = 0
        self.linecol = canvas_1.create_line( 0,0,1,0, dash = (2,2) )
        self.text_distint = canvas_1.create_text(0,0,text="",font=("arial",8))
    def integrate(self):
        self.o1_point = eval("self.o1." + str(self.o1_coo))
        self.coocol = ClosestPointOnLineSegment(self.o1_point[0],self.o1_point[1],self.o2)
        self.distint = sqrt( (self.coocol[0]-self.o1_point[0])**2 + (self.coocol[1]-self.o1_point[1])**2 )
        canvas_1.coords(self.linecol,self.coocol[0],self.coocol[1],self.o1_point[0],self.o1_point[1])
        canvas_1.coords(self.text_distint, self.coocol[0] + 0.5*(self.o1_point[0] - self.coocol[0]) , self.coocol[1] + 0.5*(self.o1_point[1] - self.coocol[1]) )
        #canvas_1.itemconfigure(self.text_distint,text=round(self.distint))

Dp = DebugPoint(0,0)

obj.append(Object_Ball(canvas_1,500,200,10,rho_rubber,0.5,0.6)) #spring cube ball 1
obj[0].v0[0] = 0.00
obj[0].v0[1] = -0.0

obj.append(Object_Ball(canvas_1,450,200,10,rho_rubber,0.5,0.6)) #spring cube ball 2
obj[1].v0[0] = 0.0
obj[1].v0[1] = 0.0

obj.append(Object_Ball(canvas_1,500,150,10,rho_rubber,0.5,0.6)) #spring cube ball 3
obj[2].v0[0] = 0.0
obj[2].v0[1] = 0.0

obj.append(Object_Ball(canvas_1,450,150,10,rho_rubber,0.5,0.6)) #spring cube ball 4
obj[3].v0[0] = 0.0
obj[3].v0[1] = 0.0

obj.append(Object_Ball(canvas_1,480,230,15,rho_rubber,0.5,0.6))
obj[4].v0[0] = 3.0
obj[4].v0[1] = 0.0

obj.append(Object_Ball(canvas_1,130,50,15,rho_steel,0.5,0.6))
obj[5].v0[0] = 0.0
obj[5].v0[1] = 0.0

obj.append(Object_Ball(canvas_1,100,15,8,rho_rubber,0.5,0.6))
obj[6].v0[0] = 0.0
obj[6].v0[1] = 0.0

obj.append(Object_Ball(canvas_1,120,55,10,rho_rubber,0.5,0.6))
obj[7].v0[0] = 0.0
obj[7].v0[1] = 0.0

obj.append(Object_FixedPoint(canvas_1,200,230))

obj.append(Object_Ball(canvas_1,140,70,12,rho_rubber,0.5,0.6))
obj[9].v0[0] = 0.0
obj[9].v0[1] = 0.0

obj.append(Object_FixedLine(canvas_1,320,200,445,300,5,0.5,0.6))
obj.append(Object_FixedLine(canvas_1,450,500,600,500,5,0.5,0.6))

obj.append(Object_FixedLine(canvas_1,315,h_cnv-100,800,h_cnv-100,5,0.5,0.3))

obj.append(Object_Trace(canvas_1,obj[5],10,20))

#mouse
obj.append(Object_FixedPoint(canvas_1,400,400))
obj.append(Object_Ball(canvas_1,400,450,14,rho_rubber,0.5,0.6))

obj.append(Object_Box(canvas_1,600,350,70,30,50,rho_rubber))

# spring cube
con.append(Constraint_SpringDamper(obj[0],obj[1],100,1000000,500000))
con.append(Constraint_SpringDamper(obj[1],obj[2],100,1000000,500000))
con.append(Constraint_SpringDamper(obj[2],obj[3],100,1000000,500000))
con.append(Constraint_SpringDamper(obj[3],obj[0],100,1000000,500000))
con.append(Constraint_SpringDamper(obj[2],obj[0],100,1000000,500000))
con.append(Constraint_SpringDamper(obj[1],obj[3],100,1000000,500000))
#

# 4 balls to 1 fixed point
con.append(Constraint_SpringDamper(obj[5],obj[8],40,1000000,1000000))
con.append(Constraint_SpringDamper(obj[8],obj[6],25,1000000,1000000))
con.append(Constraint_SpringDamper(obj[6],obj[7],125,1000000,1000000))
con.append(Constraint_SpringDamper(obj[9],obj[6],125,1000000,1000000))

# To mouse
#con.append(Constraint_SpringDamper(obj[14],obj[15],0,10000000,1000000))
con.append(Constraint_Distance(obj[14],obj[15],50))

# Rope
obj.append(Object_FixedPoint(canvas_1,650,50))
obj.append(Object_Ball(canvas_1,650,100,10,rho_rubber,0.5,0.6))
obj.append(Object_Ball(canvas_1,645,150,10,rho_rubber,0.5,0.6))
obj.append(Object_Ball(canvas_1,655,200,10,rho_rubber,0.5,0.6))
obj.append(Object_Ball(canvas_1,675,250,10,rho_rubber,0.5,0.6))
obj.append(Object_Ball(canvas_1,695,300,10,rho_rubber,0.5,0.6))
con.append(Constraint_Distance(obj[17],obj[18],50))
con.append(Constraint_Distance(obj[18],obj[19],50))
con.append(Constraint_Distance(obj[19],obj[20],50))
con.append(Constraint_Distance(obj[20],obj[21],50))
con.append(Constraint_Distance(obj[21],obj[22],50))

#obj.append(Object_ShowDistance_Point_LineSegment(obj[4],obj[10]))
#obj.append(Object_ShowDistance_Point_LineSegment(obj[4],obj[11]))

obj.append(Object_Line(canvas_1,260,150,380,140,10,40,rho_steel,0.25,0.2)) #270,150 | 390,70

obj.append(Object_ShowDistance_Point_LineSegment(obj[15],"coo1",obj[10]))
obj.append(Object_ShowDistance_Point_LineSegment(obj[15],"coo1",obj[11]))

obj.append(Object_ShowPhysics(canvas_1,obj[4]))
obj.append(Object_Trace(canvas_1,obj[4],10,50))

obj.append(Object_Ball(canvas_1,300,100,20,rho_steel,0.6,0.6)) #big ball 1
obj.append(Object_Ball(canvas_1,300,50,10,rho_steel,0.6,0.6)) #big ball 2

obj.append(Object_ShowPhysics(canvas_1,obj[28]))
obj.append(Object_ShowPhysics(canvas_1,obj[29]))

obj.append(Object_FixedLine(canvas_1,320,210,205,230,5,0.5,0.6))

#obj.append(Object_FixedLine(canvas_1,400,180,145,180,0.5,0.6))

obj.append(Object_ShowDistance_Point_LineSegment(obj[23],"cooA",obj[12]))
obj.append(Object_ShowDistance_Point_LineSegment(obj[23],"cooB",obj[12]))
obj.append(Object_ShowPhysics(canvas_1,obj[23]))

obj.append(Object_Line(canvas_1,350,125,410,115,10,60,rho_rubber,0.5,0.5)) #270,150 | 390,70


x_m = 0
y_m = 0
def motion(event):
    global x_m, y_m
    x_m, y_m = event.x, event.y
    #print('{}, {}'.format(x_m, y_m))
    obj[14].coofix = [x_m, y_m]
    canvas_1.coords(obj[14].ball,x_m-obj[14].r,y_m-obj[14].r,x_m+obj[14].r,y_m+obj[14].r)
canvas_1.bind('<Motion>', motion)

#Testkommentar för githubcommit

## ATT GÖRA:
# Varför funkar den nya impulsekvationen i L2FL och inte den gamla?
# Ha kraftpilar osv som ett separat objekt. Som Trace
# fixa initial velocity för ball
# cube 
# för fyrkanter/avlånga objekt osv, beräkna dess tvärsnittsarea som möter vinden varje tidssteg - påverkar luftmotsståndet - kan skapa vingar osv
# skapa musdrag av prylar
# skapa rep, med inställning av antalet delar och längd
# allmän kodförbättring - förbättra integrationen så att den är mer samlad?
# skapa portaler?


Timestep()

win_1.mainloop()
