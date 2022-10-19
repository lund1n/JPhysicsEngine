#!/usr/bin/python
from tkinter import *
from numpy import *
import time

from numpy import arcsin

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
rho_honey = 1400.0
rho_steel = 8000.0
rho_rubber = 920.0
rho_medium = rho_air

obj = []
con = []
viz_sf = 40
viz_sf2 = 1000000
sf_farrow1 = 50
sf_farrow2 = 100000000
viz = []

colcounter = 0

debug_col_lines = 0

t_elapsed = 0
dt = 0.035 #0.002 0.015 0.05

dtr = canvas_1.create_text(10,20,text="dt = "+str(dt),font=("arial",8),anchor=SW)
t_elapsedr = canvas_1.create_text(10,30,text="t = "+str(round(t_elapsed,9)),font=("arial",8),anchor=SW)

linecol = canvas_1.create_line( 0,0,1,0, dash = (2,2) ) # vad är detta?

# coordinate system
# x
canvas_1.create_line(20,50,60,50,arrow=LAST,fill="red")
canvas_1.create_text(66,50,text="x",font=("arial",8),fill="red")
# y
canvas_1.create_line(20,50,20,90,arrow=LAST,fill="green")
canvas_1.create_text(20,96,text="y",font=("arial",8),fill="green")

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
            #print(obj[i].F_D_form)
        #if i == 15: # the mouse ball
        #    print(obj[i].F_D_form)
    
    #print(obj[4].coo0[1]-obj[4].coom1[1])
    #print( sqrt( (obj[4].coo0[0]-obj[4].coom1[0])**2+(obj[4].coo0[1]-obj[4].coom1[1])**2 ) )

    # Loop
    win_1.after(5,Timestep)

def Rotvec2(v,ang):
    return array([ cos(ang)*v[0] - sin(ang)*v[1] , sin(ang)*v[0] + cos(ang)*v[1] ])

def Rotxypairsinvec(v,ang):
    output = []
    for i in range(0,len(v),2):
        output.append(cos(ang)*v[i] - sin(ang)*v[i+1])
        output.append(sin(ang)*v[i] + cos(ang)*v[i+1])
    return output

def Rotvec2_90degcw(v):
    return array([ v[1] , -v[0] ])

def Rotvec2_90degcountercw(v):
    return array([ -v[1] , v[0] ])

def Rotzvec3_90degcw(v):
    return array([ v[1] , -v[0] , v[2] ])

def Safediv(up,down):
    if down != 0:
        return up/down
    #except TypeError:
    #    return [ up[0]/down , up[1]/down ]
    else:
        return sign(up)*100000000

def Normalize(vec2):
    return Safediv( vec2 , sqrt( vec2[0]**2 + vec2[1]**2 ) )

def Magnitude(vec2):
    return sqrt( vec2[0]**2 + vec2[1]**2 )

def Crossprod_vec2(v,w):
    return v[0]*w[1]-v[1]*w[0]

def Crossprod_vec3(v,w):
    return [ v[1]*w[2] - v[2]*w[1] , v[2]*w[0] - v[0]*w[2] , v[0]*w[1] - v[1]*w[0] ]

def Matmultcol(M,v):
    return [M[0][0]*v[0] , M[1][1]*v[1] , M[2][2]*v[2]]

def Dotprod_vec2(v,w):
    return v[0]*w[0]+v[1]*w[1]

def Scalarproj_vec3(v,w):
    return Safediv( (v[0]*w[0] + v[1]*w[1] + v[2]*w[2]) , sqrt(w[0]**2 + w[1]**2 + w[2]**2) )

def Vecproj_vec2(v,w):
    return multiply( Safediv( v[0]*w[0] + v[1]*w[1] , (w[0]**2 + w[1]**2) ) , w )

def Vecproj_vec3(v,w):
    return multiply( Safediv( v[0]*w[0] + v[1]*w[1] + v[2]*w[2] , (w[0]**2 + w[1]**2 + w[2]**2) ) , w )

### ej längre nödvändig
def ClosestPointOnLineSegment(px,py,line):
    '''
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
    '''

    theta_line = line.theta1
    offset = [ line.coovertex[0] , line.coovertex[1] ]

    # remove offset from origin
    p = subtract([px,py],offset)
    vector = subtract( [ line.coovertex[2] , line.coovertex[3] ],offset)

    # transform from global coordinates x and y, to local coordinates u and v
    p = Rotvec2(p,-theta_line)
    vector = Rotvec2(vector,-theta_line)

    vecumax = max(vector[0],0)
    vecumin = min(vector[0],0)

    # project vertically
    p[1] = 0

    # check edge cases
    if p[0] < vecumin:
        p = [vecumin,0]
    if p[0] > vecumax:
        p = [vecumax,0]

    # transform back
    p = Rotvec2(p,theta_line)
    # add offset from origin
    p = add(p,offset)

    #return cooint
    return p
###

def ClosestPointOnLineSegmentAbsDist(px,py,line): #WORK HERE - mata in angle så du slipper räkna ut det varje gång?
    '''
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
    '''
    theta_obj = arctan( Safediv( line[3]-line[1] , line[2]-line[0] ) )
    #theta_obj = obj.theta1
    offset = [ line[0] , line[1] ]

    # remove offset from origin
    p = subtract([px,py],offset)
    vector = subtract( [ line[2] , line[3] ],offset)

    # transform from global coordinates x and y, to local coordinates u and v
    p = Rotvec2(p,-theta_obj)
    vector = Rotvec2(vector,-theta_obj)

    distint = abs(p[1])

    vecumax = max(vector[0],0)
    vecumin = min(vector[0],0)

    # p_outside_of_line_edges = point projection perpendicular to vector, outside of vector bounds
    p_outside_of_line_edges = 0

    # check edge cases
    if p[0] < vecumin:
        distint = sqrt( (p[0]-vecumin)**2 + (p[1])**2 )
        p = [vecumin,0]
        p_outside_of_line_edges = 1
    if p[0] > vecumax:
        distint = sqrt( (p[0]-vecumax)**2 + (p[1])**2 )
        p = [vecumax,0]
        p_outside_of_line_edges = 1

    # project vertically
    p[1] = 0

    # transform back
    p = Rotvec2(p,theta_obj)
    # add offset from origin
    p = add(p,offset)

    #return cooint
    return [ p[0] , p[1] , distint, p_outside_of_line_edges ]

def ClosestPointOnLineSegmentPerpDist(px,py,line,gcx,gcy,t): # WORK HERE
    '''
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
    '''


    theta_obj = arctan( Safediv( line[3]-line[1] , line[2]-line[0] ) )
    #theta_obj = obj.theta1
    offset = [ line[0] , line[1] ]

    # remove offset from origin
    p = subtract([px,py],offset)
    vector = subtract( [ line[2] , line[3] ],offset)
    gc = subtract([gcx,gcy],offset)

    # transform from global coordinates x and y, to local coordinates u and v
    p = Rotvec2(p,-theta_obj)
    vector = Rotvec2(vector,-theta_obj)
    gc = Rotvec2(gc,-theta_obj)

    #angle = arctan(Safediv( (line[3]-line[1]) , (line[2]-line[0]) ))*180/pi

    distpen = p[1]*sign(gc[1])+t # >0 = inside, <0 = outside
    #if distpen > 0 and (abs(p[1])>abs(gc[1])):
    #    distpen = -1

    dist_dyngc_int = abs(p[1])

    #sign_distint_p_signed = p[1]
    #sign_int_gc_signed = gc[1]
    #p_gc_opposing_sides_of_line = 0
    #if sign_distint_p_signed != sign_int_gc_signed:
    #    p_gc_opposing_sides_of_line = abs(gc[1])
    
    vecumax = max(vector[0],0)
    vecumin = min(vector[0],0)

    # p_outside_of_line_edges = point projection perpendicular to vector, outside of vector bounds
    p_outside_of_line_edges = 0

    # check edge cases
    if p[0] < vecumin:
        #distint = sqrt( (p[0]-vecumin)**2 + (p[1])**2 )
        #p = [vecumin,0]
        #distpen = -2
        p_outside_of_line_edges = 1
    if p[0] > vecumax:
        #distint = sqrt( (p[0]-vecumax)**2 + (p[1])**2 )
        #p = [vecumax,0]
        #distpen = -3
        p_outside_of_line_edges = 1
    if (sign(p[1]) == sign(gc[1])) and (abs(p[1])>abs(gc[1])):
        #distpen = -4
        p_outside_of_line_edges = 1
    

    kA = gc[1]/(gc[0]-vecumin)
    mA = gc[1]-kA*gc[0]
    xlimA = (p[1]-mA)/kA
    kB = (-gc[1])/(vecumax-gc[0])
    mB = gc[1]-kB*gc[0]
    xlimB = (p[1]-mB)/kB
    #xlimB = p[1]/kB+gc[0]
    #if p_outside_of_line_edges == 0:
    
    #unitvec_n = [0,-sign(gc[1])]
    '''
    if distpen>0:
        ofsx = 100
        ofsy = 500
        #viz.append( canvas_1.create_rectangle(ofs-100,ofs-100,ofs+100,ofs+100,fill="white") )
        ( canvas_1.create_line(vecumin+ofsx,ofsy,100+ofsx,ofsy,fill="black",arrow=LAST) ) #x axis
        ( canvas_1.create_line(vecumin+ofsx,ofsy,vecumin+ofsx,100+ofsy,fill="black",arrow=LAST) ) #y axis
        
        ( canvas_1.create_line(vecumin+ofsx,ofsy,vecumax+ofsx,ofsy,fill="red") )
        ( canvas_1.create_line(vecumax+ofsx,ofsy,vecumax+ofsx,2*gc[1]+ofsy,fill="red") )
        ( canvas_1.create_line(vecumin+ofsx,ofsy,vecumin+ofsx,2*gc[1]+ofsy,fill="red") )
        ( canvas_1.create_line(vecumin+ofsx,2*gc[1]+ofsy,vecumax+ofsx,2*gc[1]+ofsy,fill="red") )

        ( canvas_1.create_line(gc[0]+ofsx,gc[1]+ofsy,gc[0]+ofsx,gc[1]+ofsy+unitvec_n[1]*25,fill="grey",arrow=LAST) ) # unitvec_n

        ( canvas_1.create_oval(p[0]-3+ofsx,p[1]-3+ofsy,p[0]+3+ofsx,p[1]+3+ofsy,fill="purple",outline="purple") )
        #viz.append( canvas_1.create_line(p[0]+ofsx,ofsy,p[0]+ofsx,p[1]+ofsy,fill="purple") )
        #viz.append( canvas_1.create_oval(p[0]-3+ofsx,-3+ofsy,p[0]+3+ofsx,3+ofsy,fill="red",outline="red") )

        ( canvas_1.create_oval(gc[0]-3+ofsx,gc[1]-3+ofsy,gc[0]+3+ofsx,gc[1]+3+ofsy,fill="orange",outline="orange") )
        ( canvas_1.create_line(vecumin+ofsx,ofsy,gc[0]+ofsx,gc[1]+ofsy,fill="orange") )
        ( canvas_1.create_line(gc[0]+ofsx,gc[1]+ofsy,vecumax+ofsx,ofsy,fill="orange") )
    '''
    #print(distpen)
    
    #if distpen > 0:
    #    #print(distpen)
    if distpen > 0 and (p[0]) < xlimA:
        #print("px = "+str(p[0])+" < xlimA = "+str(xlimA))
        #print("py = "+str(p[1]))
        p_outside_of_line_edges = 1
    elif distpen > 0 and (p[0]) > xlimB:
        #print("px = "+str(p[0])+" > xlimB = "+str(xlimB))
        #print("py = "+str(p[1]))
        p_outside_of_line_edges = 1
    #else:
    #    print("px = "+str(p[0])+", py = "+str(p[1])+", xlimA = "+str(xlimA)+", xlimB = "+str(xlimB))
        #print("#################")
    
    #unitvec_n = Rotvec2(unitvec_n,theta_obj)
    #unitvec_n = [unitvec_n[0],unitvec_n[1],0]

    # project vertically
    p[1] = 0 #- sign(gc[1])*t
    
    # transform back
    p = Rotvec2(p,theta_obj)
    # add offset from origin
    p = add(p,offset)

    #if p_outside_of_line_edges == 0:
    #    if distpen > 0:
    #        viz.append( canvas_1.create_line(p[0],p[1],px,py,fill="green") )
    #    if distpen == -4:
    #        viz.append( canvas_1.create_line(p[0],p[1],px,py,fill="orange",dash = (2,2)) )
    #else:
    #    viz.append( canvas_1.create_line(p[0],p[1],px,py,fill="red") )

    #return cooint

    #testdist = sqrt( (px-p[0])**2 + (py-p[1])**2)

    return [ p[0] , p[1] , distpen, p_outside_of_line_edges, dist_dyngc_int ]


def ClosestPointOnLineSegmentPerpDist2(px,py,line,angle): # WORK HERE
    
    theta_obj = angle
    #theta_obj = arctan( Safediv( line[3]-line[1] , line[2]-line[0] ) )
    offset = [ line[0] , line[1] ]

    # remove offset from origin
    p = subtract([px,py],offset)
    vector = subtract( [ line[2] , line[3] ],offset)

    # transform from global coordinates x and y, to local coordinates u and v
    p = Rotvec2(p,-theta_obj)
    vector = Rotvec2(vector,-theta_obj)

    dist_to_surface = abs(p[1])

    vecumax = max(vector[0],0)
    vecumin = min(vector[0],0)

    p_outside_of_line_edges = 0

    # check edge cases
    if p[0] < vecumin:
        p_outside_of_line_edges = 1
    if p[0] > vecumax:
        p_outside_of_line_edges = 1
        
    # project vertically
    p[1] = 0
    
    # transform back
    p = Rotvec2(p,theta_obj)
    # add offset from origin
    p = add(p,offset)

    return [ p[0] , p[1] , dist_to_surface, p_outside_of_line_edges ]

def Pointinsidepolygoncheck(pointcoo,polygoncoo):

    #inside_or_outside_checklist = (len(polygoncoo)/2)*[0] # create a list entry for every line in the polygon
    point_inside_polygon = 0
    ints_h = 0
    ints_v = 0

    for i in range(int(len(polygoncoo)/2)):
        if i == (int(len(polygoncoo)/2)-1): # if i is at the last xy-pair in the list
            line = [polygoncoo[2*i],polygoncoo[2*i+1],polygoncoo[0],polygoncoo[1]]
        else: # else, count as usual
            line = [polygoncoo[2*i],polygoncoo[2*i+1],polygoncoo[2*i+2],polygoncoo[2*i+3]]

        line_max_x = max(line[0],line[2])
        line_min_x = min(line[0],line[2])
        line_max_y = max(line[1],line[3])
        line_min_y = min(line[1],line[3])

        if line[1] == line[3]: # is the line horizontal?
            # raycast upwards
            if pointcoo[0] >= line_min_x and pointcoo[0] < line_max_x and line[1] <= pointcoo[1]:
                ints_v = ints_v + 1

        elif line[0] == line[2]: # is the line vertical?
            # raycast to the sides
            if pointcoo[1] >= line_min_y and pointcoo[1] <= line_max_y and line[0] <= pointcoo[0]:
                ints_h = ints_h + 1

        else: # now we can safely check the line which must be sloped
            k = (line[3]-line[1]) / (line[2]-line[0]) # already determined to not be zero or infinity. Correct
            m = line[1]-k*line[0]

            #x_int_h = (pointcoo[1]-m)/k
            #y_int_h = k*x_int_h+m
            
            #y_int_v = k*pointcoo[0]+m
            #x_int_v = (y_int_v-m)/k

            y_int_v = k*pointcoo[0]+m
            x_int_v = (y_int_v-m)/k

            ##################################################################
            ### DEBUG DRAW SECTION ###
            ##################################################################
            if pointcoo[0] == obj[40].coovertex[4] and i==1:
                viz.append( canvas_1.create_oval(x_int_v-3,y_int_v-3,x_int_v+3,y_int_v+3,fill="orange",outline="black") )
                viz.append( canvas_1.create_line( obj[40].coovertex[4],0,obj[40].coovertex[4],800, dash = (2,2) ) )
                
            if pointcoo[0] == obj[40].coovertex[4] and i==3:
                viz.append( canvas_1.create_oval(x_int_v-3,y_int_v-3,x_int_v+3,y_int_v+3,fill="purple",outline="black") )
            ##################################################################

            #canvas_1.coords(p2_draw,x_int_h-3,y_int_h-3,x_int_h+3,y_int_h+3)
            #canvas_1.coords(p3_draw,x_int_v-3,y_int_v-3,x_int_v+3,y_int_v+3)
            
            # raycast upwards
            if x_int_v >= line_min_x and x_int_v <= line_max_x and y_int_v <= pointcoo[1]:
                ints_v = ints_v + 1
            # raycast to the sides
            #if x_int_h >= line_min_x and x_int_h <= line_max_x and x_int_h <= pointcoo[0]:
            #    ints_h = ints_h + 1
            
    #if (ints_h%2) != 0 and (ints_v%2) != 0:
    #if (ints_h%2) != 0:
    if (ints_v%2) != 0:
        point_inside_polygon = 1
        #canvas_1.itemconfigure(p_draw,fill="green",outline="green")
        #canvas_1.coords(p_draw,x_m-5,y_m-5,x_m+5,y_m+5)
        #canvas_1.create_oval(x_m-2,y_m-2,x_m+2,y_m+2,fill="black",outline="black")
        
    return point_inside_polygon

def Ballinsidepolygoncheck(polygon,pointcoo,polygoncoo,ballradius,linenormals): #UNUSED AND DOES NOT WORK WELL DUE TO CONCAVE POLYGONS ETC.

    point_inside_polygon = 1

    for i in range(polygon.n_linesegments):
        if i == (polygon.n_linesegments-1): # if i is at the last xy-pair in the list
            line = [polygoncoo[2*i],polygoncoo[2*i+1],polygoncoo[0],polygoncoo[1]]
        else: # else, count as usual
            line = [polygoncoo[2*i],polygoncoo[2*i+1],polygoncoo[2*i+2],polygoncoo[2*i+3]]

        point = [ pointcoo[0]-line[0]+ballradius*linenormals[2*i] , pointcoo[1]-line[1]+ballradius*linenormals[2*i+1] ]
        dp = Dotprod_vec2( point , [linenormals[2*i],linenormals[2*i+1]] )

        if dp < 0:
            point_inside_polygon = 0

    return point_inside_polygon

def BallinsidepolygoncheckOLD(polygon,pointcoo,polygoncoo,ballradius,lineangles,linenormals):

    #inside_or_outside_checklist = (len(polygoncoo)/2)*[0] # create a list entry for every line in the polygon
    point_inside_polygon = 0
    ints_h = 0
    ints_v = 0

    for i in range(polygon.n_linesegments):
        if i == (polygon.n_linesegments-1): # if i is at the last xy-pair in the list
            line = [polygoncoo[2*i],polygoncoo[2*i+1],polygoncoo[0],polygoncoo[1]]
        else: # else, count as usual
            line = [polygoncoo[2*i],polygoncoo[2*i+1],polygoncoo[2*i+2],polygoncoo[2*i+3]]

        line_max_x = max(line[0],line[2])
        line_min_x = min(line[0],line[2])
        line_max_y = max(line[1],line[3])
        line_min_y = min(line[1],line[3])

        if line[1] == line[3]: # is the line horizontal?
            # raycast upwards
            if pointcoo[0] >= line_min_x and pointcoo[0] < line_max_x and (pointcoo[1]+ballradius*linenormals[2*i+1]) > line[1]:
                ints_v = ints_v + 1

        elif line[0] == line[2]: # is the line vertical?
            # raycast to the sides
            if pointcoo[1] >= line_min_y and pointcoo[1] < line_max_y and (pointcoo[0]+ballradius*linenormals[2*i]) > line[0]:
                ints_h = ints_h + 1

        else: # now we can safely check the line which must be sloped
            k = (line[3]-line[1]) / (line[2]-line[0]) # already determined to not be zero or infinity
            m = line[1]-k*line[0]

            x_int_h = (pointcoo[1]-m)/k
            y_int_h = k*x_int_h+m
            
            y_int_v = k*pointcoo[0]+m
            x_int_v = (y_int_v-m)/k

            #canvas_1.coords(p2_draw,x_int_h-3,y_int_h-3,x_int_h+3,y_int_h+3)
            #canvas_1.coords(p3_draw,x_int_v-3,y_int_v-3,x_int_v+3,y_int_v+3)
            # raycast upwards
            if (x_int_v-ballradius*linenormals[2*i]) >= line_min_x and (x_int_v-ballradius*linenormals[2*i]) < line_max_x and (y_int_v-ballradius*linenormals[2*i+1]) < pointcoo[1]:
                ints_v = ints_v + 1
                #print(ballradius)
                #print("#######")
                #print(y_int_v)
                #print(ballradius*abs(cos(lineangles[i])))
                #print(y_int_v+ballradius*abs(cos(lineangles[i])))
                #print(pointcoo[1])
                #print("#######")
            # raycast to the sides
            if (x_int_h-ballradius*linenormals[2*i]) >= line_min_x and (x_int_h-ballradius*linenormals[2*i]) < line_max_x and (x_int_h-ballradius*linenormals[2*i]) < pointcoo[0]:
                ints_h = ints_h + 1
                #print(ballradius*abs(sin(lineangle[i])))

            
    if (ints_h%2) != 0 and (ints_v%2) != 0:
        point_inside_polygon = 1
        #canvas_1.itemconfigure(p_draw,fill="green",outline="green")
        #canvas_1.coords(p_draw,x_m-5,y_m-5,x_m+5,y_m+5)
        #canvas_1.create_oval(x_m-2,y_m-2,x_m+2,y_m+2,fill="black",outline="black")
        
    return point_inside_polygon

'''
#OLD CODE
def IntersectionPoint_Line_Line(l1,l2):
    cooint = [0,0]
    # 2 line intersection
    cooint[0] = Safediv( ( (l2.cooA[1]-l2.k*l2.cooA[0])-(l1.cooA[1]-l1.k*l1.cooA[0]) ),(l1.k-l2.k) ) #xi = (m2-m1)/(k1-k2) = ( (y2-k2x2)-(y1-k1x1) )/(k1-k2)
    cooint[1] = l1.k*cooint[0]+(l1.cooA[1]-l1.k*l1.cooA[0]) # yi = k1xi+m1 = k1xi+(y1-k1x1)
    return cooint
'''
def IntersectionPoint_Line_Line(l1,l2):

    iols = 0 # intersection inside of line segment

    if l1[0] == l1[2]: # line 1 vertical
        if l2[0] == l2[2]: # line 2 vertical
            return [ 0 , 0 , 0 ]
        elif l2[1] == l2[3]: # line 2 horizontal
            iols = ( l1[0]>=min(l2[0],l2[2]) )*( l1[0]<=max(l2[0],l2[2]) )*( min(l1[1],l1[3])<=l2[1] )*( max(l1[1],l1[3])>=l2[1] )
            return [ l1[0] , l2[1] , iols ]
        else: # line 2 angled
            k2 = (l2[3]-l2[1]) / (l2[2]-l2[0]) # already determined to not be zero or infinity
            m2 = l2[1]-k2*l2[0]
            x_int = l1[0]
            y_int = k2*x_int+m2
            iols = ( x_int>=min(l2[0],l2[2]) )*( x_int<=max(l2[0],l2[2]) )*( min(l1[1],l1[3])<=y_int )*( max(l1[1],l1[3])>=y_int )
            return [ x_int , y_int , iols ]

    elif l1[1] == l1[3]: # line 1 horizontal
        if l2[0] == l2[2]: # line 2 vertical
            iols = ( l2[0]>=min(l1[0],l1[2]) )*( l2[0]<=max(l1[0],l1[2]) )*( min(l2[1],l2[3])<=l1[1] )*( max(l2[1],l2[3])>=l1[1] )
            return [ l2[0] , l1[1] , iols ]
        elif l2[1] == l2[3]: # line 2 horizontal
            return [ 0 , 0 , 0 ]
        else: # line 2 angled
            k2 = (l2[3]-l2[1]) / (l2[2]-l2[0]) # already determined to not be zero or infinity
            m2 = l2[1]-k2*l2[0]
            y_int = l1[1]
            x_int = (y_int-m2)/k2
            iols = ( x_int>=min(l1[0],l1[2]) )*( x_int<=max(l1[0],l1[2]) )*( min(l2[1],l2[3])<=y_int )*( max(l2[1],l2[3])>=y_int )
            return [ x_int , y_int , iols ]
    else: # line 1 angled
        k1 = (l1[3]-l1[1]) / (l1[2]-l1[0]) # already determined to not be zero or infinity
        m1 = l1[1]-k1*l1[0]
        if l2[0] == l2[2]: # line 2 vertical
            x_int = l2[0]
            y_int = k1*x_int+m1
            iols = ( x_int>=min(l1[0],l1[2]) )*( x_int<=max(l1[0],l1[2]) )*( min(l2[1],l2[3])<=y_int )*( max(l2[1],l2[3])>=y_int )
            return [ x_int , y_int , iols ]
        elif l2[1] == l2[3]: # line 2 horizontal
            y_int = l2[1]
            x_int = (y_int-m1)/k1
            iols = ( x_int>=min(l2[0],l2[2]) )*( x_int<=max(l2[0],l2[2]) )*( min(l1[1],l1[3])<=y_int )*( max(l1[1],l1[3])>=y_int )
            return [ x_int , y_int , iols ]
        else:
            k2 = (l2[3]-l2[1]) / (l2[2]-l2[0]) # already determined to not be zero or infinity
            m2 = l2[1]-k2*l2[0]
            x_int = (m2-m1)/(k1-k2)
            y_int = k1*x_int+m1
            iols = ( x_int>=min(l2[0],l2[2]) )*( x_int<=max(l2[0],l2[2]) )*( min(l1[1],l1[3])<=y_int )*( max(l1[1],l1[3])>=y_int )
            return [ x_int , y_int , iols ]

def Check_CollisionType(o1,o2):

    if isinstance(o1, Object_Ball):
        if isinstance(o2, Object_Ball):
            Collision_Ball_Ball(o1,o2)
            return
        if isinstance(o2, Object_Polygon):
            Collision_Ball_Polygon(o1,o2)
            return
        if isinstance(o2, Object_Line):
            Collision_Ball_Line(o1,o2)
            return
        #if obj[4] is o1:
        #if 1 == 1:
        if isinstance(o2, Object_FixedLine):
            Collision_Ball_FixedLine(o1,o2)
            return
        #if isinstance(o2, Object_Polygon):
        #    Collision_Ball_Box(o2,o1)
        #    return

    if isinstance(o1, Object_FixedLine):
        if isinstance(o2, Object_Ball):
            Collision_Ball_FixedLine(o2,o1)
            return
        if isinstance(o2, Object_Line):
            Collision_Line_FixedLine(o2,o1)
            return
        if isinstance(o2, Object_Polygon):
            Collision_Polygon_FixedLine(o2,o1)
            return

    if isinstance(o1, Object_Line):
        if isinstance(o2, Object_FixedLine):
            Collision_Line_FixedLine(o1,o2)
            return
        if isinstance(o2, Object_Ball):
            Collision_Ball_Line(o2,o1)
            return
        if isinstance(o2, Object_Line):
            Collision_Line_Line(o2,o1)
            return
        if isinstance(o2, Object_Polygon):
            Collision_Box_Line(o2,o1)
            return

    if isinstance(o1, Object_Polygon):
        if isinstance(o2, Object_FixedLine):
            Collision_Polygon_FixedLine(o1,o2)
            return
        if isinstance(o2, Object_Line):
            Collision_Box_Line(o1,o2)
            return
        #if isinstance(o2, Object_Ball):
        #    Collision_Ball_Box(o1,o2)
        #    return
        if isinstance(o2, Object_Ball):
            Collision_Ball_Polygon(o2,o1)
            return
        if isinstance(o2, Object_Polygon):
            Collision_Polygon_Polygon(o1,o2)
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


    # Draw normal arrows
    #viz.append( canvas_1.create_line(coocol[0]-unitvec_n[0]*viz_sf,coocol[1]-unitvec_n[1]*viz_sf,coocol[0]+unitvec_n[0]*viz_sf,coocol[1]+unitvec_n[1]*viz_sf,arrow=BOTH,fill="green") )
    viz.append( canvas_1.create_line(coocol[0]-unitvec_n[0]*sf_farrow1*tanh(abs(F_n)/sf_farrow2),coocol[1]-unitvec_n[1]*sf_farrow1*tanh(abs(F_n)/sf_farrow2),coocol[0]+unitvec_n[0]*sf_farrow1*tanh(abs(F_n)/sf_farrow2),coocol[1]+unitvec_n[1]*sf_farrow1*tanh(abs(F_n)/sf_farrow2),arrow=BOTH,fill="red") )
    # Draw tang arrows
    #viz.append( canvas_1.create_line(coocol[0]-unitvec_t[0]*30,coocol[1]-unitvec_t[1]*30,coocol[0]+unitvec_t[0]*30,coocol[1]+unitvec_t[1]*30,arrow=BOTH,fill="green") )
    viz.append( canvas_1.create_line(coocol[0]-unitvec_t[0]*sf_farrow1*tanh(abs(F_t)/sf_farrow2),coocol[1]-unitvec_t[1]*sf_farrow1*tanh(abs(F_t)/sf_farrow2),coocol[0]+unitvec_t[0]*sf_farrow1*tanh(abs(F_t)/sf_farrow2),coocol[1]+unitvec_t[1]*sf_farrow1*tanh(abs(F_t)/sf_farrow2),arrow=BOTH,fill="orange") )
    # Draw contact point
    viz.append( canvas_1.create_oval(coocol[0]-3,coocol[1]-3,coocol[0]+3,coocol[1]+3,fill="red",outline="red") )

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
    viz.append( canvas_1.create_oval(coocol[0]-3,coocol[1]-3,coocol[0]+3,coocol[1]+3,fill="red",outline="red") )
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
    viz.append( canvas_1.create_oval(coocol[0]-3,coocol[1]-3,coocol[0]+3,coocol[1]+3,fill="red",outline="red") )
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
        b2.coom1 = subtract( b2.coom1 , multiply(b2.m/(b1.m+b2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))


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
    coocol = ClosestPointOnLineSegmentAbsDist(ball.coo1[0],ball.coo1[1],line.coovertex)
    # Distance between object and intersection point
    distint2 = sqrt( (coocol[0]-ball.coo1[0])**2 + (coocol[1]-ball.coo1[1])**2 ) # distance between point and fixline
    distint = coocol[2]

    if distint < ball.r+line.hhalf:
        unitvec_n = divide( array([ ball.coo1[0]-coocol[0] , ball.coo1[1]-coocol[1] , 0 ]) , distint ) # Unit normal vector
        unitvec_t = Rotzvec3_90degcw(unitvec_n) # Unit tangent vector
        #unitvec_n = line.unitvec_n
        #unitvec_t = line.unitvec_t

        coocol_hcomp = add( [coocol[0],coocol[1],0], multiply( line.hhalf, [ unitvec_n[0], unitvec_n[1] , 0 ] ) )

        Contact_dyn_statline(ball,line,coocol_hcomp,unitvec_n,unitvec_t)

        # Old. DO NOT REMOVE. 
        #ball.F_other = add( ball.F_other , multiply( nvec , (2*ball.m/(dt**2))*(ball.r-distint) ) )

        ball.coo1 = ball.coo1 + multiply([unitvec_n[0],unitvec_n[1]],(ball.r+line.hhalf-distint))
        ball.coo0 = ball.coo0 + multiply([unitvec_n[0],unitvec_n[1]],(ball.r+line.hhalf-distint))
        ball.coom1 = ball.coom1 + multiply([unitvec_n[0],unitvec_n[1]],(ball.r+line.hhalf-distint))

def Collision_Ball_Line(ball,line):
    # Coordinates of intersection point
    coocol = ClosestPointOnLineSegmentAbsDist(ball.coo1[0],ball.coo1[1],line.coovertex)
    # Distance between object and intersection point
    #distint = sqrt( (coocol[0]-ball.coo1[0])**2 + (coocol[1]-ball.coo1[1])**2 ) # distance between point and fixline
    distint = coocol[2]

    if distint < ball.r+line.hhalf:
        unitvec_n = divide( array([ ball.coo1[0]-coocol[0] , ball.coo1[1]-coocol[1] , 0 ]) , distint ) # Unit normal vector
        unitvec_t = Rotzvec3_90degcw(unitvec_n) # Unit tangent vector
    
        coocol_hcomp = add( [coocol[0],coocol[1],0], multiply( line.hhalf, [ unitvec_n[0], unitvec_n[1], 0 ] ) )

        Contact_dyn_line(ball,line,coocol_hcomp,unitvec_n,unitvec_t) #FELFELFELFELFELFEL FIXA!!!

        pendist = -(distint-(ball.r+line.hhalf))

        ball.coo1 = add( ball.coo1 , multiply(ball.m/(ball.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
        ball.coo0 = add( ball.coo0 , multiply(ball.m/(ball.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
        ball.coom1 = add( ball.coom1 , multiply(ball.m/(ball.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
        line.coo1 = subtract( line.coo1 , multiply(line.m/(ball.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
        line.coo0 = subtract( line.coo0 , multiply(line.m/(ball.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
        line.coom1 = subtract( line.coom1 , multiply(line.m/(ball.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
        
        #vert1 = subtract( [line.coovertex[0],line.coovertex[1]] , multiply(line.m/(ball.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
        #vert2 = subtract( [line.coovertex[2],line.coovertex[3]] , multiply(line.m/(ball.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
        #line.coovertex = [vert1[0],vert1[1],vert2[0],vert1[1]]

        #line.cooA = subtract( line.cooA , multiply(line.m/(ball.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
        #line.cooB = subtract( line.cooB , multiply(line.m/(ball.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))

        line.coovertex = subtract( line.coovertex , multiply(line.m/(ball.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]] ))

def Collision_Line_Line(l1,l2):
    
    r_pinball = l1.hhalf + l2.hhalf
    #fac = 1

    if sqrt( (l1.coo0[0]-l2.coo0[0])**2 + (l1.coo0[1]-l2.coo0[1])**2 ) <= (l1.r + l2.r + r_pinball):

        coocol_l1A = ClosestPointOnLineSegmentAbsDist(l1.coovertex[0],l1.coovertex[1],l2.coovertex)
        coocol_l1B = ClosestPointOnLineSegmentAbsDist(l1.coovertex[2],l1.coovertex[3],l2.coovertex)
        coocol_l2A = ClosestPointOnLineSegmentAbsDist(l2.coovertex[0],l2.coovertex[1],l1.coovertex)
        coocol_l2B = ClosestPointOnLineSegmentAbsDist(l2.coovertex[2],l2.coovertex[3],l1.coovertex)
        distint_l1A = coocol_l1A[2]
        distint_l1B = coocol_l1B[2]
        distint_l2A = coocol_l2A[2]
        distint_l2B = coocol_l2B[2]

        if debug_col_lines==1:
            viz.append( canvas_1.create_line(coocol_l1A[0],coocol_l1A[1],coocol_l2A[0],coocol_l2A[1],fill="red") )
            viz.append( canvas_1.create_line(coocol_l1A[0],coocol_l1A[1],coocol_l2B[0],coocol_l2B[1],fill="green") )
            viz.append( canvas_1.create_line(coocol_l1B[0],coocol_l1B[1],coocol_l2A[0],coocol_l2A[1],fill="blue") )
            viz.append( canvas_1.create_line(coocol_l1B[0],coocol_l1B[1],coocol_l2B[0],coocol_l2B[1],fill="purple") )

        lines_edge_to_edge = coocol_l1A[3]*coocol_l1B[3]*coocol_l2A[3]*coocol_l2B[3]

        #distint_l1A = sqrt( (coocol_l1A[0]-l1.coovertex[0])**2 + (coocol_l1A[1]-l1.coovertex[1])**2 )
        #distint_l1B = sqrt( (coocol_l1B[0]-l1.coovertex[2])**2 + (coocol_l1B[1]-l1.coovertex[3])**2 )
        #distint_l2A = sqrt( (coocol_l2A[0]-l2.coovertex[0])**2 + (coocol_l2A[1]-l2.coovertex[1])**2 )
        #distint_l2B = sqrt( (coocol_l2B[0]-l2.coovertex[2])**2 + (coocol_l2B[1]-l2.coovertex[3])**2 )
        
        #print(l2)
        #print(obj[38])
        #if l2 == obj[38]:
        #    print(1)

        #if distint_l1A <= r_pinball and (coocol_l1A[3]==0 or (coocol_l1A[3]==1 and (coocol_l2A[3]==1 or coocol_l2B[3]==1))):
        if distint_l1A <= r_pinball and (coocol_l1A[3]==0 or lines_edge_to_edge==1):
            unitvec_n = divide( [ l1.coovertex[0] - coocol_l1A[0] , l1.coovertex[1] - coocol_l1A[1] , 0 ] , distint_l1A ) # pekar alltid mot linje 1
            unitvec_t = Rotzvec3_90degcw(unitvec_n)

            coocol_hcomp = add( [coocol_l1A[0],coocol_l1A[1],0], multiply( l2.hhalf, unitvec_n ) )

            Contact_line_line(l1,l2,coocol_hcomp,unitvec_n,unitvec_t)
    
            pendist = (r_pinball-distint_l1A)
            l1.coo1 = add( l1.coo1 , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l1.coo0 = add( l1.coo0 , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l1.coom1 = add( l1.coom1 , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            #l1.cooA = add( l1.cooA , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            #l1.cooB = add( l1.cooB , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l1.coovertex = add( l1.coovertex , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]] ))

            l2.coo1 = subtract( l2.coo1 , multiply(l2.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l2.coo0 = subtract( l2.coo0 , multiply(l2.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l2.coom1 = subtract( l2.coom1 , multiply(l2.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l2.coovertex = subtract( l2.coovertex , multiply(l2.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]] ))

        #if distint_l1B <= r_pinball and (coocol_l1B[3]==0 or (coocol_l1B[3]==1 and (coocol_l2A[3]==1 or coocol_l2B[3]==1))):
        if distint_l1B <= r_pinball and (coocol_l1B[3]==0 or lines_edge_to_edge==1):
            unitvec_n = divide( [ l1.coovertex[2] - coocol_l1B[0] , l1.coovertex[3] - coocol_l1B[1] , 0 ] , distint_l1B ) # pekar alltid mot linje 1
            unitvec_t = Rotzvec3_90degcw(unitvec_n)

            coocol_hcomp = add( [coocol_l1B[0],coocol_l1B[1],0], multiply( l2.hhalf, unitvec_n ) )

            Contact_line_line(l1,l2,coocol_hcomp,unitvec_n,unitvec_t)
        
            pendist = (r_pinball-distint_l1B)
            l1.coo1 = add( l1.coo1 , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l1.coo0 = add( l1.coo0 , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l1.coom1 = add( l1.coom1 , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            #l1.cooA = add( l1.cooA , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            #l1.cooB = add( l1.cooB , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l1.coovertex = add( l1.coovertex , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]] ))
        
            l2.coo1 = subtract( l2.coo1 , multiply(l2.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l2.coo0 = subtract( l2.coo0 , multiply(l2.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l2.coom1 = subtract( l2.coom1 , multiply(l2.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l2.coovertex = subtract( l2.coovertex , multiply(l2.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]] ))
        
        if distint_l2A <= r_pinball and (coocol_l2A[3]==0 or lines_edge_to_edge==1):
            unitvec_n = divide( [ coocol_l2A[0] - l2.coovertex[0] , coocol_l2A[1] - l2.coovertex[1] , 0 ] , distint_l2A ) # pekar alltid mot linje 1
            unitvec_t = Rotzvec3_90degcw(unitvec_n)

            coocol_hcomp = subtract( [coocol_l2A[0],coocol_l2A[1],0], multiply( l1.hhalf, unitvec_n ) )
            
            Contact_line_line(l1,l2,coocol_hcomp,unitvec_n,unitvec_t)

            pendist = (r_pinball-distint_l2A)
            l1.coo1 = add( l1.coo1 , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l1.coo0 = add( l1.coo0 , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l1.coom1 = add( l1.coom1 , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            #l1.cooA = add( l1.cooA , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            #l1.cooB = add( l1.cooB , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l2.coovertex = add( l1.coovertex , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]] ))
            
            l2.coo1 = subtract( l2.coo1 , multiply(l2.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l2.coo0 = subtract( l2.coo0 , multiply(l2.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l2.coom1 = subtract( l2.coom1 , multiply(l2.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l2.coovertex = subtract( l2.coovertex , multiply(l2.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]] ))

        if distint_l2B <= r_pinball and (coocol_l2B[3]==0 or lines_edge_to_edge==1):
            unitvec_n = divide( [ coocol_l2B[0] - l2.coovertex[2] , coocol_l2B[1] - l2.coovertex[3] , 0 ] , distint_l2B ) # pekar alltid mot linje 1
            unitvec_t = Rotzvec3_90degcw(unitvec_n)

            coocol_hcomp = subtract( [coocol_l2B[0],coocol_l2B[1],0], multiply( l1.hhalf, unitvec_n ) )

            Contact_line_line(l1,l2,coocol_hcomp,unitvec_n,unitvec_t)
            
            pendist = (r_pinball-distint_l2B)
            l1.coo1 = add( l1.coo1 , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l1.coo0 = add( l1.coo0 , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l1.coom1 = add( l1.coom1 , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            #l1.cooA = add( l1.cooA , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            #l1.cooB = add( l1.cooB , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l1.coovertex = add( l1.coovertex , multiply(l1.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]] ))
            
            l2.coo1 = subtract( l2.coo1 , multiply(l2.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l2.coo0 = subtract( l2.coo0 , multiply(l2.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l2.coom1 = subtract( l2.coom1 , multiply(l2.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
            l2.coovertex = subtract( l2.coovertex , multiply(l2.m/(l1.m+l2.m)*pendist, [unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]] ))

def Collision_Line_FixedLine(line,fl):
    
    r_pinball = line.hhalf + fl.hhalf
    #fac = 1

    #global colcounter

    if sqrt( (line.coo0[0]-fl.coo0[0])**2 + (line.coo0[1]-fl.coo0[1])**2 ) <= (line.r + fl.r + r_pinball):

        coocol_box2_lineA = ClosestPointOnLineSegmentAbsDist(line.coovertex[0],line.coovertex[1],fl.coovertex)
        coocol_box2_lineB = ClosestPointOnLineSegmentAbsDist(line.coovertex[2],line.coovertex[3],fl.coovertex)
        coocol_flA = ClosestPointOnLineSegmentAbsDist(fl.coovertex[0],fl.coovertex[1],line.coovertex)
        coocol_flB = ClosestPointOnLineSegmentAbsDist(fl.coovertex[2],fl.coovertex[3],line.coovertex)
        distint_lineA = coocol_box2_lineA[2]
        distint_lineB = coocol_box2_lineB[2]
        distint_flA = coocol_flA[2]
        distint_flB = coocol_flB[2]
        #distint_lineA = sqrt( (coocol_box2_lineA[0]-line.coovertex[0])**2 + (coocol_box2_lineA[1]-line.coovertex[1])**2 )
        #distint_lineB = sqrt( (coocol_box2_lineB[0]-line.coovertex[2])**2 + (coocol_box2_lineB[1]-line.coovertex[3])**2 )
        #distint_flA = sqrt( (coocol_flA[0]-fl.coovertex[0])**2 + (coocol_flA[1]-fl.coovertex[1])**2 )
        #distint_flB = sqrt( (coocol_flB[0]-fl.coovertex[2])**2 + (coocol_flB[1]-fl.coovertex[3])**2 )

        if debug_col_lines==1:
            viz.append( canvas_1.create_line(coocol_box2_lineA[0],coocol_box2_lineA[1],line.coovertex[0],line.coovertex[1],fill="green") )
            viz.append( canvas_1.create_line(coocol_box2_lineB[0],coocol_box2_lineB[1],line.coovertex[2],line.coovertex[3],fill="green") )
            viz.append( canvas_1.create_line(coocol_flA[0],coocol_flA[1],fl.coovertex[0],fl.coovertex[1],fill="green") )
            viz.append( canvas_1.create_line(coocol_flB[0],coocol_flB[1],fl.coovertex[2],fl.coovertex[3],fill="green") )

        # a line can only collide with another line at two locations
        #cols_detected = 0
        lines_edge_to_edge = coocol_box2_lineA[3]*coocol_box2_lineB[3]*coocol_flA[3]*coocol_flB[3]

        #if distint_lineA <= r_pinball and (coocol_box2_lineA[3]==0 or (coocol_box2_lineA[3]==1 and (coocol_flA[3]==1 or coocol_flB[3]==1))):
        
        if distint_lineA <= r_pinball and (coocol_box2_lineA[3]==0 or lines_edge_to_edge==1):

            #if line == obj[36]:
            #    colcounter = colcounter + 1
            #    print(str(colcounter)+" - lineA - "+str(coocol_box2_lineA[2])+" ( "+str(coocol_box2_lineA[0])+" , "+str(coocol_box2_lineA[1])+" )")

            unitvec_n = divide( [ line.coovertex[0] - coocol_box2_lineA[0] , line.coovertex[1] - coocol_box2_lineA[1] , 0 ] , distint_lineA )
            unitvec_t = Rotzvec3_90degcw(unitvec_n)

            coocol_hcomp = add( [coocol_box2_lineA[0],coocol_box2_lineA[1],0], multiply( fl.hhalf, unitvec_n ) )

            Contact_dyn_statline(line,fl,coocol_hcomp,unitvec_n,unitvec_t)
    
            line.coo1 = line.coo1 + multiply([unitvec_n[0],unitvec_n[1]],(r_pinball-distint_lineA))
            line.coo0 = line.coo0 + multiply([unitvec_n[0],unitvec_n[1]],(r_pinball-distint_lineA))
            line.coom1 = line.coom1 + multiply([unitvec_n[0],unitvec_n[1]],(r_pinball-distint_lineA))
            #line.cooA = line.cooA + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_lineA))
            #line.cooB = line.cooB + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_lineA))
            line.coovertex = line.coovertex + multiply([unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]],(r_pinball-distint_lineA))

            #cols_detected = cols_detected + 1
    
        #if distint_lineB <= r_pinball and (coocol_box2_lineB[3]==0 or (coocol_box2_lineB[3]==1 and (coocol_flA[3]==1 or coocol_flB[3]==1))):
        if distint_lineB <= r_pinball and (coocol_box2_lineB[3]==0 or lines_edge_to_edge==1):

            #if line == obj[36]:
            #    colcounter = colcounter + 1
            #    print(str(colcounter)+" - lineB - "+str(coocol_box2_lineB[2])+" ( "+str(coocol_box2_lineB[0])+" , "+str(coocol_box2_lineB[1])+" )")

            unitvec_n = divide( [ line.coovertex[2] - coocol_box2_lineB[0] , line.coovertex[3] - coocol_box2_lineB[1] , 0 ] , distint_lineB )
            unitvec_t = Rotzvec3_90degcw(unitvec_n)

            coocol_hcomp = add( [coocol_box2_lineB[0],coocol_box2_lineB[1],0], multiply( fl.hhalf, unitvec_n ) )

            Contact_dyn_statline(line,fl,coocol_hcomp,unitvec_n,unitvec_t)
        
            line.coo1 = line.coo1 + multiply([unitvec_n[0],unitvec_n[1]],(r_pinball-distint_lineB))
            line.coo0 = line.coo0 + multiply([unitvec_n[0],unitvec_n[1]],(r_pinball-distint_lineB))
            line.coom1 = line.coom1 + multiply([unitvec_n[0],unitvec_n[1]],(r_pinball-distint_lineB))
            #line.cooA = line.cooA + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_lineB))
            #line.cooB = line.cooB + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_lineB))
            line.coovertex = line.coovertex + multiply([unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]],(r_pinball-distint_lineB))
            
            #cols_detected = cols_detected + 1
        
        if distint_flA <= r_pinball and (coocol_flA[3]==0 or lines_edge_to_edge==1):

            #if line == obj[36]:
            #    colcounter = colcounter + 1
            #    print(str(colcounter)+" - flA - "+str(coocol_flA[2])+" ( "+str(coocol_flA[0])+" , "+str(coocol_flA[1])+" )")

            unitvec_n = divide( [ coocol_flA[0] - fl.coovertex[0] , coocol_flA[1] - fl.coovertex[1] , 0 ] , distint_flA )
            unitvec_t = Rotzvec3_90degcw(unitvec_n)

            coocol_hcomp = subtract( [coocol_flA[0],coocol_flA[1],0], multiply( line.hhalf, unitvec_n ) )
            
            Contact_dyn_statline(line,fl,coocol_hcomp,unitvec_n,unitvec_t)

            line.coo1 = line.coo1 + multiply([unitvec_n[0],unitvec_n[1]],(r_pinball-distint_flA))
            line.coo0 = line.coo0 + multiply([unitvec_n[0],unitvec_n[1]],(r_pinball-distint_flA))
            line.coom1 = line.coom1 + multiply([unitvec_n[0],unitvec_n[1]],(r_pinball-distint_flA))
            #line.cooA = line.cooA + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_flA))
            #line.cooB = line.cooB + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_flA))
            line.coovertex = line.coovertex + multiply([unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]],(r_pinball-distint_flA))

            #cols_detected = cols_detected + 1
    
        if distint_flB <= r_pinball and (coocol_flB[3]==0 or lines_edge_to_edge==1):

            #if line == obj[36]:
            #    colcounter = colcounter + 1
            #    print(str(colcounter)+" - flB - "+str(coocol_flB[2])+" ( "+str(coocol_flB[0])+" , "+str(coocol_flB[1])+" )")

            unitvec_n = divide( [ coocol_flB[0] - fl.coovertex[2] , coocol_flB[1] - fl.coovertex[3] , 0 ] , distint_flB )
            unitvec_t = Rotzvec3_90degcw(unitvec_n)

            coocol_hcomp = subtract( [coocol_flB[0],coocol_flB[1],0], multiply( line.hhalf, unitvec_n ) )

            Contact_dyn_statline(line,fl,coocol_hcomp,unitvec_n,unitvec_t)
            
            line.coo1 = line.coo1 + multiply([unitvec_n[0],unitvec_n[1]],(r_pinball-distint_flB))
            line.coo0 = line.coo0 + multiply([unitvec_n[0],unitvec_n[1]],(r_pinball-distint_flB))
            line.coom1 = line.coom1 + multiply([unitvec_n[0],unitvec_n[1]],(r_pinball-distint_flB))
            #line.cooA = line.cooA + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_flB))
            #line.cooB = line.cooB + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_flB))
            line.coovertex = line.coovertex + multiply([unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]],(r_pinball-distint_flB))

            #cols_detected = cols_detected + 1
                
        
        #if line == obj[36]:
        #    print("#############################################")
            
        #if obj[10] is fl:
        #    #print(coocol_box2_lineA)
        #    Dp.moveto( coocol_box2_lineA[0],coocol_box2_lineA[1] )


def Collision_Ball_Polygon(ball,polygon):

    if sqrt( (ball.coo0[0]-polygon.coo0[0])**2 + (ball.coo0[1]-polygon.coo0[1])**2 ) <= (ball.r + polygon.r):

        len_polygon = 2*polygon.n_linesegments
        len_polygon_half = polygon.n_linesegments

        #bip = Ballinsidepolygoncheck(polygon,[ball.coo1[0],ball.coo1[1]],polygon.coovertex,ball.r,polygon.angles_local_rotated,polygon.normals_local_rotated)
        #the bip below seems to NOT work
        #bip = Ballinsidepolygoncheck(polygon,[ball.coo1[0],ball.coo1[1]],polygon.coovertex,ball.r,polygon.normals_local_rotated)
        #the bip below works
        bip = Pointinsidepolygoncheck([ball.coo1[0],ball.coo1[1]],polygon.coovertex)
        
        polygon_unitvec_n = polygon.normals_local_rotated
        polygon_unitvec_t = polygon.tangents_local_rotated

        #if bip > 0: # if a point/vertex of ball is inside of polygon
        if 1==1:

            #print(1)
            #print(2)
            #print("---")

            point = [ball.coo1[0],ball.coo1[1]]

            dist_list = [None]*len_polygon_half
            int_list = [None]*len_polygon
            ool_list = [None]*len_polygon_half

            #####
            for j in range(len_polygon_half): # ... go through all lines in polygon
                if j == (len_polygon_half-1): # if i is at the last xy-pair in the list
                    line = [polygon.coovertex[2*j],polygon.coovertex[2*j+1],polygon.coovertex[0],polygon.coovertex[1]]
                else: # else, count as usual
                    line = [polygon.coovertex[2*j],polygon.coovertex[2*j+1],polygon.coovertex[2*j+2],polygon.coovertex[2*j+3]]
                
                line_angle = polygon.angles_local_rotated[j]

                cooint = ClosestPointOnLineSegmentAbsDist(point[0],point[1],line)
                #cooint = ClosestPointOnLineSegmentPerpDist2(point[0],point[1],line,line_angle)
                dist_list[j] = cooint[2] # store distance of closest point of o1 point on the o2 line
                int_list[2*j] = cooint[0] # store x coordinate of closest point of o1 point on the o2 line
                int_list[2*j+1] = cooint[1] # store y coordinate of closest point of o1 point on the o2 line
                ool_list[j] = cooint[3] # outside of line list. checks if the o1 point is outside of the endpoints of the o2 line and thus cant hit it
            #####

            shortest_dist = min(dist_list)
            shortest_dist_index = dist_list.index(shortest_dist)
            colpx = int_list[2*shortest_dist_index]
            colpy = int_list[2*shortest_dist_index+1]
            #colpx = colpx + polygon_unitvec_n[2*shortest_dist_index]*ball.r
            #colpy = colpy + polygon_unitvec_n[2*shortest_dist_index+1]*ball.r
            
            '''
            print("dist_list")
            print(dist_list)
            print("shortest dist")
            print(shortest_dist)
            print("shortest dist index")
            print(shortest_dist_index)
            print("ool_list")
            print(ool_list)
            print("ool_list shortest distance index")
            print(ool_list[shortest_dist_index])
            print("---------------------")
            '''

            #if ool_list[shortest_dist_index]==0:
            #if 1==1:
            if shortest_dist <= ball.r or bip:
                print("shortest dist")
                print(shortest_dist)
                print("ball.r")
                print(ball.r)
                print("ool_list shortest distance index")
                print(ool_list[shortest_dist_index])
                #canvas_1.create_oval(colpx-3,colpy-3,colpx+3,colpy+3,fill="green",outline="green")
                Contact_line_line(ball,polygon,[colpx,colpy,0],[polygon_unitvec_n[2*shortest_dist_index],polygon_unitvec_n[2*shortest_dist_index+1],0],[polygon_unitvec_t[2*shortest_dist_index],polygon_unitvec_t[2*shortest_dist_index+1],0]) # WORK HERE

                '''
                print("###############################")
                #print("i: "+str(i))
                print(shortest_dist_index)
                #print(colp_o1)
                print(dist_list)
                print(shortest_dist)
                print([unitvec_n[2*shortest_dist_index],unitvec_n[2*shortest_dist_index+1],0])
                print([unitvec_t[2*shortest_dist_index],unitvec_t[2*shortest_dist_index+1],0])
                print("colp_ball:   "+str(colp_ball))
                print("colp_polygon:   "+str(colp_polygon))
                #if o1 == ball:
                #    print("ball")
                #if o1 == polygon:
                #    print("polygon")
                print("###############################")
                '''
                
                pendist = ball.r - shortest_dist
                ball_xy_translate = multiply(ball.m/(ball.m+polygon.m)*pendist,[polygon_unitvec_n[2*shortest_dist_index],polygon_unitvec_n[2*shortest_dist_index+1]])
                polygon_xy_translate = multiply(polygon.m/(ball.m+polygon.m)*pendist,[polygon_unitvec_n[2*shortest_dist_index],polygon_unitvec_n[2*shortest_dist_index+1]])
                ball.coo1 = subtract( ball.coo1 , ball_xy_translate )
                ball.coo0 = subtract( ball.coo0 , ball_xy_translate )
                ball.coom1 = subtract( ball.coom1 , ball_xy_translate )
                
                polygon.coo1 = add( polygon.coo1 , polygon_xy_translate )
                polygon.coo0 = add( polygon.coo0 , polygon_xy_translate )
                polygon.coom1 = add( polygon.coom1 , polygon_xy_translate )
                polygon.coo_geom_center = add( polygon.coo_geom_center , polygon_xy_translate )
                for i in range(len_polygon_half):
                    polygon.coovertex[2*i] = add( polygon.coovertex[2*i] , polygon_xy_translate[0] )
                    polygon.coovertex[2*i+1] = add( polygon.coovertex[2*i+1] , polygon_xy_translate[1] )
                        
def Collision_Ball_Box(box,ball):
    
    #r_pinball = ball.r

    if sqrt( (box.coo0[0]-ball.coo0[0])**2 + (box.coo0[1]-ball.coo0[1])**2 ) <= (box.r + ball.r):

        indices = [ [0,1,2,3] , [2,3,4,5] , [4,5,6,7] , [6,7,0,1] ]

        for i in range(0,len(indices)):
            boxline_active = [box.coovertex[ indices[i][0] ],box.coovertex[ indices[i][1] ],box.coovertex[ indices[i][2] ],box.coovertex[ indices[i][3] ]]

            #coocol_ball = ClosestPointOnLineSegmentAbsDist(ball.coo1[0],ball.coo1[1],boxline_active)
            #distint_ball = coocol_ball[2]

            coocol_ball = ClosestPointOnLineSegmentPerpDist(ball.coo1[0],ball.coo1[1],boxline_active,box.coo_geom_center[0],box.coo_geom_center[1],ball.r)
            distpen_ball = coocol_ball[2]
            dist_ballgc_int = coocol_ball[4]

            #if distpen_ball <= r_pinball and coocol_ball[3]==0:
            if (distpen_ball) > 0 and coocol_ball[3]==0:
                #unitvec_n = divide( [ coocol_ball[0] - ball.coo1[0] , coocol_ball[1] - ball.coo1[1] , 0 ] , dist_ballgc_int )
                unitvec_n = divide( [ ball.coo1[0] - coocol_ball[0] , ball.coo1[1] - coocol_ball[1] , 0 ] , dist_ballgc_int )
                #unitvec_n = divide( [ coocol_ball[0] - ball.coo1[0] , coocol_ball[1] - ball.coo1[1] , 0 ] , sqrt( (coocol_ball[0] - ball.coo1[0])**2 + (coocol_ball[1] - ball.coo1[1])**2) )

                #unitvec_n = divide( [ ball.coo1[0] - coocol_ball[0] , ball.coo1[1] - coocol_ball[1] , 0 ] , distpen_ball )
                #unitvec_n = divide( array([ ball.coo1[0]-coocol_ball[0] , ball.coo1[1]-coocol_ball[1] , 0 ]) , distint_ball ) # Unit normal vector
                unitvec_t = Rotzvec3_90degcw(unitvec_n)
                #print(unitvec_n)
                #print( sqrt( (unitvec_n[0])**2 + (unitvec_n[1])**2 ) )
                #print( sqrt( (coocol_ball[0] - ball.coo1[0])**2 + (coocol_ball[1] - ball.coo1[1])**2 ) )
                #print([ coocol_ball[0] - ball.coo1[0] , coocol_ball[1] - ball.coo1[1] , 0 ])
                #print(distpen_ball)
                
                #canvas_1.create_oval(coocol_ball[0]-3,coocol_ball[1]-3,coocol_ball[0]+3,coocol_ball[1]+3,fill="blue",outline="blue")

                Contact_dyn_line(ball,box,coocol_ball,unitvec_n,unitvec_t)

                pendist = distpen_ball
                box.coo1 = subtract( box.coo1 , multiply(box.m/(box.m+ball.m)*pendist,[unitvec_n[0],unitvec_n[1]]) )
                box.coo0 = subtract( box.coo0 , multiply(box.m/(box.m+ball.m)*pendist,[unitvec_n[0],unitvec_n[1]]) )
                box.coom1 = subtract( box.coom1 , multiply(box.m/(box.m+ball.m)*pendist,[unitvec_n[0],unitvec_n[1]]) )
                #box.coovertex = subtract( box.coovertex , multiply(box.m/(box.m+ball.m)*pendist,[unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]]) )
                for j in range(box.n_linesegments):
                    box.coovertex[2*j] = box.coovertex[2*j] - unitvec_n[0]*(box.m/(box.m+ball.m)*pendist)
                    box.coovertex[2*j+1] = box.coovertex[2*j+1] - unitvec_n[1]*(box.m/(box.m+ball.m)*pendist)

                ball.coo1 = add( ball.coo1 , multiply(ball.m/(box.m+ball.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
                ball.coo0 = add( ball.coo0 , multiply(ball.m/(box.m+ball.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
                ball.coom1 = add( ball.coom1 , multiply(ball.m/(box.m+ball.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))

def Collision_Ball_Box_OLD(box,ball):
    
    r_pinball = ball.r

    if sqrt( (box.coo0[0]-ball.coo0[0])**2 + (box.coo0[1]-ball.coo0[1])**2 ) <= (box.r + ball.r):

        indices = [ [0,1,2,3] , [2,3,4,5] , [4,5,6,7] , [6,7,0,1] ]

        for i in range(0,len(indices)):
            boxline_active = [box.coovertex[ indices[i][0] ],box.coovertex[ indices[i][1] ],box.coovertex[ indices[i][2] ],box.coovertex[ indices[i][3] ]]

            coocol_ball = ClosestPointOnLineSegmentAbsDist(ball.coo1[0],ball.coo1[1],boxline_active)
            distint_ball = coocol_ball[2]

            if distint_ball <= r_pinball and coocol_ball[3]==0:
                unitvec_n = divide( array([ ball.coo1[0]-coocol_ball[0] , ball.coo1[1]-coocol_ball[1] , 0 ]) , distint_ball ) # Unit normal vector
                #unitvec_n = divide( [ boxline_active[0] - coocol_ball[0] , boxline_active[1] - coocol_ball[1] , 0 ] , distint_ball )
                #unitvec_n = ( [ Safediv(coocol_ball[0] - box.coovertex[0] , distint_ball) , Safediv(coocol_ball[1] - box.coovertex[1] , distint_ball) , 0 ] )
                unitvec_t = Rotzvec3_90degcw(unitvec_n)

                Contact_dyn_line(ball,box,coocol_ball,unitvec_n,unitvec_t)

                pendist = (distint_ball-r_pinball)
                box.coo1 = add( box.coo1 , multiply(box.m/(box.m+ball.m)*pendist,[unitvec_n[0],unitvec_n[1]]) )
                box.coo0 = add( box.coo0 , multiply(box.m/(box.m+ball.m)*pendist,[unitvec_n[0],unitvec_n[1]]) )
                box.coom1 = add( box.coom1 , multiply(box.m/(box.m+ball.m)*pendist,[unitvec_n[0],unitvec_n[1]]) )
                box.coovertex = add( box.coovertex , multiply(box.m/(box.m+ball.m)*pendist,[unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]]) )

                ball.coo1 = subtract( ball.coo1 , multiply(ball.m/(box.m+ball.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
                ball.coo0 = subtract( ball.coo0 , multiply(ball.m/(box.m+ball.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
                ball.coom1 = subtract( ball.coom1 , multiply(ball.m/(box.m+ball.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))

def Collision_Box_Line(box,line):
    
    r_pinball = line.hhalf

    if sqrt( (box.coo0[0]-line.coo0[0])**2 + (box.coo0[1]-line.coo0[1])**2 ) <= (box.r + line.r + r_pinball):

        #indices = [ [0,1,2,3] , [2,3,4,5] , [4,5,6,7] , [6,7,0,1] ]

        for i in range(0,box.n_linesegments):
            #boxline_active = [box.coovertex[ indices[i][0] ],box.coovertex[ indices[i][1] ],box.coovertex[ indices[i][2] ],box.coovertex[ indices[i][3] ]]
            if i == (box.n_linesegments-1): # if i is at the last xy-pair in the list
                boxline_active = [box.coovertex[2*i],box.coovertex[2*i+1],box.coovertex[0],box.coovertex[1]]
            else: # else, count as usual
                boxline_active = [box.coovertex[2*i],box.coovertex[2*i+1],box.coovertex[2*i+2],box.coovertex[2*i+3]]

            coocol_lineA = ClosestPointOnLineSegmentAbsDist(line.coovertex[0],line.coovertex[1],boxline_active)
            coocol_lineB = ClosestPointOnLineSegmentAbsDist(line.coovertex[2],line.coovertex[3],boxline_active)
            distint_lineA = coocol_lineA[2]
            distint_lineB = coocol_lineB[2]

            coocol_boxA = ClosestPointOnLineSegmentAbsDist(boxline_active[0],boxline_active[1],line.coovertex)
            coocol_boxB = ClosestPointOnLineSegmentAbsDist(boxline_active[2],boxline_active[3],line.coovertex)
            distint_boxA = coocol_boxA[2]
            distint_boxB = coocol_boxB[2]

            if debug_col_lines==1:
                viz.append( canvas_1.create_line(coocol_boxA[0],coocol_boxA[1],boxline_active[0],boxline_active[1],fill="green") )
                viz.append( canvas_1.create_line(coocol_boxB[0],coocol_boxB[1],boxline_active[2],boxline_active[3],fill="green") )
                viz.append( canvas_1.create_line(coocol_lineA[0],coocol_lineA[1],line.coovertex[0],line.coovertex[1],fill="green") )
                viz.append( canvas_1.create_line(coocol_lineB[0],coocol_lineB[1],line.coovertex[2],line.coovertex[3],fill="green") )

            if distint_boxA <= r_pinball and coocol_boxA[3]==0:
                unitvec_n = divide( [ boxline_active[0] - coocol_boxA[0] , boxline_active[1] - coocol_boxA[1] , 0 ] , distint_boxA )
                unitvec_t = Rotzvec3_90degcw(unitvec_n)

                coocol_hcomp = add( [coocol_boxA[0],coocol_boxA[1],0], multiply( line.hhalf, unitvec_n ) )

                Contact_dyn_line(box,line,coocol_hcomp,unitvec_n,unitvec_t)
        
                pendist = (r_pinball-distint_boxA)
                box.coo1 = add( box.coo1 , multiply(box.m/(box.m+line.m)*pendist,[unitvec_n[0],unitvec_n[1]]) )
                box.coo0 = add( box.coo0 , multiply(box.m/(box.m+line.m)*pendist,[unitvec_n[0],unitvec_n[1]]) )
                box.coom1 = add( box.coom1 , multiply(box.m/(box.m+line.m)*pendist,[unitvec_n[0],unitvec_n[1]]) )
                box.coo_geom_center = add( box.coo_geom_center , multiply(box.m/(box.m+line.m)*pendist,[unitvec_n[0],unitvec_n[1]]) )
                #box.coovertex = add( box.coovertex , multiply(box.m/(box.m+line.m)*pendist,[unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]]) )
                for j in range(box.n_linesegments):
                    box.coovertex[2*j] = box.coovertex[2*j] + unitvec_n[0]*(box.m/(box.m+line.m)*pendist)
                    box.coovertex[2*j+1] = box.coovertex[2*j+1] + unitvec_n[1]*(box.m/(box.m+line.m)*pendist)

                line.coo1 = subtract( line.coo1 , multiply(line.m/(box.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
                line.coo0 = subtract( line.coo0 , multiply(line.m/(box.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
                line.coom1 = subtract( line.coom1 , multiply(line.m/(box.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
                line.coovertex = subtract( line.coovertex , multiply(line.m/(box.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]] ))
                

            if distint_lineA <= r_pinball and coocol_lineA[3]==0:
                unitvec_n = divide( [ coocol_lineA[0] - line.coovertex[0] , coocol_lineA[1] - line.coovertex[1] , 0 ] , distint_lineA )
                unitvec_t = Rotzvec3_90degcw(unitvec_n)

                Contact_dyn_line(box,line,coocol_lineA,unitvec_n,unitvec_t)

                pendist = (r_pinball-distint_lineA)
                box.coo1 = add( box.coo1 , multiply(box.m/(box.m+line.m)*pendist,[unitvec_n[0],unitvec_n[1]]) )
                box.coo0 = add( box.coo0 , multiply(box.m/(box.m+line.m)*pendist,[unitvec_n[0],unitvec_n[1]]) )
                box.coom1 = add( box.coom1 , multiply(box.m/(box.m+line.m)*pendist,[unitvec_n[0],unitvec_n[1]]) )
                box.coo_geom_center = add( box.coo_geom_center , multiply(box.m/(box.m+line.m)*pendist,[unitvec_n[0],unitvec_n[1]]) )
                #box.coovertex = add( box.coovertex , multiply(box.m/(box.m+line.m)*pendist,[unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]]) )
                for j in range(box.n_linesegments):
                    box.coovertex[2*j] = box.coovertex[2*j] + unitvec_n[0]*(box.m/(box.m+line.m)*pendist)
                    box.coovertex[2*j+1] = box.coovertex[2*j+1] + unitvec_n[1]*(box.m/(box.m+line.m)*pendist)

                line.coo1 = subtract( line.coo1 , multiply(line.m/(box.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
                line.coo0 = subtract( line.coo0 , multiply(line.m/(box.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
                line.coom1 = subtract( line.coom1 , multiply(line.m/(box.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
                line.coovertex = subtract( line.coovertex , multiply(line.m/(box.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]] ))

            if distint_lineB <= r_pinball and coocol_lineB[3]==0:
                unitvec_n = divide( [ coocol_lineB[0] - line.coovertex[2] , coocol_lineB[1] - line.coovertex[3] , 0 ] , distint_lineB )
                unitvec_t = Rotzvec3_90degcw(unitvec_n)

                Contact_dyn_line(box,line,coocol_lineB,unitvec_n,unitvec_t)

                pendist = (r_pinball-distint_lineB)
                box.coo1 = add( box.coo1 , multiply(box.m/(box.m+line.m)*pendist,[unitvec_n[0],unitvec_n[1]]) )
                box.coo0 = add( box.coo0 , multiply(box.m/(box.m+line.m)*pendist,[unitvec_n[0],unitvec_n[1]]) )
                box.coom1 = add( box.coom1 , multiply(box.m/(box.m+line.m)*pendist,[unitvec_n[0],unitvec_n[1]]) )
                box.coo_geom_center = add( box.coo_geom_center , multiply(box.m/(box.m+line.m)*pendist,[unitvec_n[0],unitvec_n[1]]) )
                #box.coovertex = add( box.coovertex , multiply(box.m/(box.m+line.m)*pendist,[unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]]) )
                for j in range(box.n_linesegments):
                    box.coovertex[2*j] = box.coovertex[2*j] + unitvec_n[0]*(box.m/(box.m+line.m)*pendist)
                    box.coovertex[2*j+1] = box.coovertex[2*j+1] + unitvec_n[1]*(box.m/(box.m+line.m)*pendist)

                line.coo1 = subtract( line.coo1 , multiply(line.m/(box.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
                line.coo0 = subtract( line.coo0 , multiply(line.m/(box.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
                line.coom1 = subtract( line.coom1 , multiply(line.m/(box.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
                line.coovertex = subtract( line.coovertex , multiply(line.m/(box.m+line.m)*pendist, [unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]] ))

def Collision_Polygon_Polygon(pg1,pg2):

    if sqrt( (pg1.coo0[0]-pg2.coo0[0])**2 + (pg1.coo0[1]-pg2.coo0[1])**2 ) <= (pg1.r + pg2.r): #if inside bounding box (change this to not use sqrt)
        intp_pg1 = 0
        intp_pg2 = 0
        #len_pg1 = len(pg1.coovertex)
        #len_pg1_half = int(len(pg1.coovertex)*0.5) #changed from /2 to *0.5. may or may not cause issues?
        #len_pg2 = len(pg2.coovertex)
        #len_pg2_half = int(len(pg2.coovertex)*0.5) #changed from /2 to *0.5. may or may not cause issues?

        len_pg1 = 2*pg1.n_linesegments
        len_pg1_half = pg1.n_linesegments
        len_pg2 = 2*pg2.n_linesegments
        len_pg2_half = pg2.n_linesegments

        colp_pg1 = [0]*len_pg1_half
        colp_pg2 = [0]*len_pg2_half

        #print("angle_pg1:   "+str(pg1.theta1)+",   angle_pg2:   "+str(pg2.theta1))

        # WORK HERE PIPC NEDAN REGISTRERAR IBLAND INTE EN GENOMBRYTNING AV EN PUNKT
        # MISSAR DEN EN LINJE? SAKNAS FALLET "om sist i listan: 2*i, 2*i+1, 0, 1"?

        # Are any points of polygon 1 inside polygon 2?
        for i in range(len_pg1_half):
            check_pg1 = Pointinsidepolygoncheck([pg1.coovertex[2*i],pg1.coovertex[2*i+1]],pg2.coovertex)
            intp_pg1 = intp_pg1 + check_pg1 #if the point is intersecting, increase this counter
            colp_pg1[i] = check_pg1 #mark the vertex which is intersecting with a 1

        # Are any points of polygon 2 inside polygon 1?
        for i in range(len_pg2_half):
            check_pg2 = Pointinsidepolygoncheck([pg2.coovertex[2*i],pg2.coovertex[2*i+1]],pg1.coovertex)
            intp_pg2 = intp_pg2 + check_pg2
            colp_pg2[i] = check_pg2

        pg1_unitvec_n = pg1.normals_local_rotated
        pg1_unitvec_t = pg1.tangents_local_rotated
        pg2_unitvec_n = pg2.normals_local_rotated
        pg2_unitvec_t = pg2.tangents_local_rotated

        ###########################################################################################################################
        ##### Define collision function #####
        ###########################################################################################################################

        def Contact_o1_o2(o1,o2,n_int,lo1,lo1h,lo2,lo2h,unitvec_n,unitvec_t,colp_o1):

            if n_int > 0: # if a point/vertex of o1 is inside of o2

                ### testsection
                coolinecross_list = [None]*lo2
                linesegcrossed_o1_list = [None]*lo1h
                linesegcrossed_o2_list = [None]*lo2h

                for i in range(lo1h): # for every point in o1 ...
                    if i == (lo1h-1): # if i is at the last xy-pair in the list
                        line1 = [o1.coovertex[2*i],o1.coovertex[2*i+1],o1.coovertex[0],o1.coovertex[1]]
                    else: # else, count as usual
                        line1 = [o1.coovertex[2*i],o1.coovertex[2*i+1],o1.coovertex[2*i+2],o1.coovertex[2*i+3]]

                    for j in range(lo2h): # ... go through all lines in o2
                        if j == (lo2h-1): # if i is at the last xy-pair in the list
                            line2 = [o2.coovertex[2*j],o2.coovertex[2*j+1],o2.coovertex[0],o2.coovertex[1]]
                        else: # else, count as usual
                            line2 = [o2.coovertex[2*j],o2.coovertex[2*j+1],o2.coovertex[2*j+2],o2.coovertex[2*j+3]]

                        crosspoint = IntersectionPoint_Line_Line([line1[0],line1[1],line1[2],line1[3]],[line2[0],line2[1],line2[2],line2[3]])
                        if crosspoint[2] == 1: #if the crossing point of line1 and line2 is inside of the line2 endpoints
                            coolinecross_list[2*j] = crosspoint[0]
                            coolinecross_list[2*j+1] = crosspoint[1]
                            linesegcrossed_o1_list[i] = 1 #DEBUG ONLY (not anymore?)
                            linesegcrossed_o2_list[j] = 1 #then mark line2 as a crossed segment of polygon2
                ### end testsection

                for i in range(lo1h): # for every point in o1 ...
                    '''
                    if i == (lo1h-1): # if i is at the last xy-pair in the list
                        line1 = [o1.coovertex[2*i],o1.coovertex[2*i+1],o1.coovertex[0],o1.coovertex[1]]
                    else: # else, count as usual
                        line1 = [o1.coovertex[2*i],o1.coovertex[2*i+1],o1.coovertex[2*i+2],o1.coovertex[2*i+3]]
                    '''
                    ###########################################################################################################################
                    ##### Info necessary for the collision event #####
                    ###########################################################################################################################

                    point = [o1.coovertex[2*i],o1.coovertex[2*i+1]]

                    dist_list = [None]*lo2h
                    int_list = [None]*lo2
                    ool_list = [None]*lo2h
                    distlinecross_list = [-1]*lo2h


                    crossed = 0  #crossed amount of unique line segments (NOT amount of crossings detected)

                    #####
                    for j in range(lo2h): # ... go through all lines in o2
                        if j == (lo2h-1): # if i is at the last xy-pair in the list
                            line2 = [o2.coovertex[2*j],o2.coovertex[2*j+1],o2.coovertex[0],o2.coovertex[1]]
                        else: # else, count as usual
                            line2 = [o2.coovertex[2*j],o2.coovertex[2*j+1],o2.coovertex[2*j+2],o2.coovertex[2*j+3]]
                        
                        line2_angle = o2.angles_local_rotated[j]

                        cooint = ClosestPointOnLineSegmentPerpDist2(point[0],point[1],line2,line2_angle)
                        int_list[2*j] = cooint[0] # store x coordinate of closest point to o1 point on the o2 line
                        int_list[2*j+1] = cooint[1] # store y coordinate of closest point to o1 point on the o2 line
                        dist_list[j] = cooint[2] # store distance of closest point to o1 point on the o2 line
                        ool_list[j] = cooint[3] # outside of line list. checks if the o1 point is outside of the endpoints of the o2 line and thus cant hit it
                        
                        '''
                        crosspoint = IntersectionPoint_Line_Line([line1[0],line1[1],line1[2],line1[3]],[line2[0],line2[1],line2[2],line2[3]])
                        if crosspoint[2] == 1: #if the crossing point of line1 and line2 is inside of the line2 endpoints
                            coolinecross_list[2*j] = crosspoint[0]
                            coolinecross_list[2*j+1] = crosspoint[1]
                            linesegcrossed_o1_list[i] = 1 #DEBUG ONLY
                            linesegcrossed_o2_list[j] = 1 #then mark line2 as a crossed segment of polygon2
                            1==1
                        '''

                    #####
                    for k in range(lo1h): #which line crossing point is the closest? it determines the line whose normal shall be used
                        #distlinecross_list[k] = 
                        if linesegcrossed_o1_list[k] == 1: #if the line segment of polygon1 has been crossed by a line segment of polygon2
                            crossed = crossed + 1 #crossed amount of unique line segments (NOT amount of crossings detected)
                    for k in range(lo2h): #which line crossing point is the closest? it determines the line whose normal shall be used
                        if linesegcrossed_o2_list[k] == 1: #if the line segment of polygon2 has been crossed by a line segment of polygon1
                            cpolsei = ClosestPointOnLineSegmentAbsDist(point[0],point[1],o2.linesegment(k))
                            distlinecross_list[k] = cpolsei[2] #record the distance to the penetrating point
                            crossed = crossed + 1 #crossed amount of unique line segments (NOT amount of crossings detected)
                        #we now check if both polygons have crossed line segments (as opposed to just o2 as originally) and increment the "crossed" variable,
                        #because if the combined number of line segments crossed is less than 3, which
                        #for some reason sometimes happens when one box slides on top of another and two edge nodes pass,
                        #one line segment crossing fails to become detected - the wrong face for separation is chosen as no other alternative is available,
                        #and objects may separate by a long distance along the wrong face normal
                    
                    filtered = []

                    for l in range(lo2h):
                        if distlinecross_list[l] > -1:
                            filtered.append(distlinecross_list[l]) #WORK HERE WORK HERE WORK HERE WORK HERE
                    if len(filtered) > 0:
                        shortest_dist = min(filtered)


                    #shortest_dist = min(dist_list)
                    #shortest_dist_index = dist_list.index(shortest_dist)
                    #colpx = int_list[2*shortest_dist_index]
                    #colpy = int_list[2*shortest_dist_index+1]
                    #####unitvec_n = o2.normals_local_rotated
                    #####unitvec_t = o2.tangents_local_rotated

                    ###########################################################################################################################
                    ##### Collision event #####                    
                    ###########################################################################################################################

                    if colp_o1[i]==1: # if the current point (point i) of object1 is inside of object2
                        #if ool_list[shortest_dist_index]==0:

                        #create line to line intersection check here to see which o2 line has been penetrated WORK HERE
                        #print("point_pg1:                 "+str(i))
                        #print("lines crossed indexes:     "+str(linesegcrossed_o2_list))
                        #print("line crossing coordinates: "+str(coolinecross_list))
                        #viz.append( canvas_1.create_oval(o1.coo0[2*i]-3,o1.coo0[2*i+1]-3,o1.coo0[2*i]+3,[2*i+1]+3,fill="red",outline="red") )

                        #if 1==1:
                        if crossed>2: # ugly fix of ignoring when sometimes only one line segment crossing per object is detected (total segments less than 3, as mentioned above)
                            
                            shortest_dist_index = distlinecross_list.index(shortest_dist)
                            #colpx = coolinecross_list[2*shortest_dist_index]
                            #colpy = coolinecross_list[2*shortest_dist_index+1]
                            #colpx = int_list[2*shortest_dist_index] #WORK HERE
                            #colpy = int_list[2*shortest_dist_index+1] #WORK HERE
                            colpx = point[0]
                            colpy = point[1]

                            #print("index, shortest distance:  "+str(shortest_dist_index))
                            #print("shortest distance:         "+str(shortest_dist))
                            #print("obj1 line segm crossed:    "+str(linesegcrossed_o1_list))
                            #print("obj2 line segm crossed:    "+str(linesegcrossed_o2_list))
                            #print("------------------------------------------------------------------------")

                            ### DEBUG - WORK HERE WORK HERE WORK HERE
                            if shortest_dist > 40:
                                1==1 # Problems occur when only one line crossing is detected per object. One of the objects should always have 2.
                            ###

                            #canvas_1.create_oval(colpx-3,colpy-3,colpx+3,colpy+3,fill="green",outline="green")
                            
                            #Definition: Contact_line_line(o1,o2,coocol,unitvec_n,unitvec_t)
                            Contact_line_line(o1,o2,[colpx,colpy,0],[unitvec_n[2*shortest_dist_index],unitvec_n[2*shortest_dist_index+1],0],[unitvec_t[2*shortest_dist_index],unitvec_t[2*shortest_dist_index+1],0]) # WORK HERE

                            '''
                            print("###############################")
                            #print("i: "+str(i))
                            print(shortest_dist_index)
                            #print(colp_o1)
                            print(dist_list)
                            print(shortest_dist)
                            print([unitvec_n[2*shortest_dist_index],unitvec_n[2*shortest_dist_index+1],0])
                            print([unitvec_t[2*shortest_dist_index],unitvec_t[2*shortest_dist_index+1],0])
                            print("colp_pg1:   "+str(colp_pg1))
                            print("colp_pg2:   "+str(colp_pg2))
                            #if o1 == pg1:
                            #    print("pg1")
                            #if o1 == pg2:
                            #    print("pg2")
                            print("###############################")
                            '''
                            
                            pendist = shortest_dist
                            o1_xy_translate = multiply(o1.m/(o1.m+o2.m)*pendist,[unitvec_n[2*shortest_dist_index],unitvec_n[2*shortest_dist_index+1]])
                            o2_xy_translate = multiply(o2.m/(o1.m+o2.m)*pendist,[unitvec_n[2*shortest_dist_index],unitvec_n[2*shortest_dist_index+1]])
                            o1.coo1 = subtract( o1.coo1 , o1_xy_translate )
                            o1.coo0 = subtract( o1.coo0 , o1_xy_translate )
                            o1.coom1 = subtract( o1.coom1 , o1_xy_translate )
                            o1.coo_geom_center = subtract( o1.coo_geom_center , o1_xy_translate )
                            for m in range(lo1h):
                                o1.coovertex[2*m] = subtract( o1.coovertex[2*m] , o1_xy_translate[0] )
                                o1.coovertex[2*m+1] = subtract( o1.coovertex[2*m+1] , o1_xy_translate[1] )
                            
                            o2.coo1 = add( o2.coo1 , o2_xy_translate )
                            o2.coo0 = add( o2.coo0 , o2_xy_translate )
                            o2.coom1 = add( o2.coom1 , o2_xy_translate )
                            o2.coo_geom_center = add( o2.coo_geom_center , o2_xy_translate )
                            for m in range(lo2h):
                                o2.coovertex[2*m] = add( o2.coovertex[2*m] , o2_xy_translate[0] )
                                o2.coovertex[2*m+1] = add( o2.coovertex[2*m+1] , o2_xy_translate[1] )
                            
                    ###########################################################################################################################

                        
                        
                        
        
        Contact_o1_o2(pg1,pg2,intp_pg1,len_pg1,len_pg1_half,len_pg2,len_pg2_half,pg2_unitvec_n,pg2_unitvec_t,colp_pg1)
        Contact_o1_o2(pg2,pg1,intp_pg2,len_pg2,len_pg2_half,len_pg1,len_pg1_half,pg1_unitvec_n,pg1_unitvec_t,colp_pg2)

def Collision_Box_Box_OLD(box1,box2):

    #r_pinball = 1
    
    if sqrt( (box1.coo0[0]-box2.coo0[0])**2 + (box1.coo0[1]-box2.coo0[1])**2 ) <= (box1.r + box2.r):

        indices = [ [0,1,2,3] , [2,3,4,5] , [4,5,6,7] , [6,7,0,1] ]
        for j in range(0,len(indices)):
            box2line_active = [box2.coovertex[ indices[j][0] ],box2.coovertex[ indices[j][1] ],box2.coovertex[ indices[j][2] ],box2.coovertex[ indices[j][3] ]]
            
            for i in range(0,len(indices)):
                box1line_active = [box1.coovertex[ indices[i][0] ],box1.coovertex[ indices[i][1] ],box1.coovertex[ indices[i][2] ],box1.coovertex[ indices[i][3] ]]
                

                coocol_box2_lineA = ClosestPointOnLineSegmentPerpDist(box2line_active[0],box2line_active[1],box1line_active,box1.coo_geom_center[0],box1.coo_geom_center[1],0)
                coocol_box2_lineB = ClosestPointOnLineSegmentPerpDist(box2line_active[2],box2line_active[3],box1line_active,box1.coo_geom_center[0],box1.coo_geom_center[1],0)
                distpen_box2_lineA = coocol_box2_lineA[2]
                distpen_box2_lineB = coocol_box2_lineB[2]

                coocol_box1_lineA = ClosestPointOnLineSegmentPerpDist(box1line_active[0],box1line_active[1],box2line_active,box2.coo_geom_center[0],box2.coo_geom_center[1],0)
                coocol_box1_lineB = ClosestPointOnLineSegmentPerpDist(box1line_active[2],box1line_active[3],box2line_active,box2.coo_geom_center[0],box2.coo_geom_center[1],0)
                distpen_box1_lineA = coocol_box1_lineA[2]
                distpen_box1_lineB = coocol_box1_lineB[2]

                if debug_col_lines==1:
                    if coocol_box1_lineA[3]==0:
                        viz.append( canvas_1.create_line(coocol_box1_lineA[0],coocol_box1_lineA[1],box1line_active[0],box1line_active[1],fill="blue") )
                    else:
                        viz.append( canvas_1.create_line(coocol_box1_lineA[0],coocol_box1_lineA[1],box1line_active[0],box1line_active[1],fill="red", dash = (3,3)) )
                    if coocol_box1_lineB[3]==0:
                        viz.append( canvas_1.create_line(coocol_box1_lineB[0],coocol_box1_lineB[1],box1line_active[2],box1line_active[3],fill="blue") )
                    else:
                        viz.append( canvas_1.create_line(coocol_box1_lineB[0],coocol_box1_lineB[1],box1line_active[2],box1line_active[3],fill="red", dash = (3,3)) )
                    if coocol_box2_lineA[3]==0:
                        viz.append( canvas_1.create_line(coocol_box2_lineA[0],coocol_box2_lineA[1],box2line_active[0],box2line_active[1],fill="blue") )
                    else:
                        viz.append( canvas_1.create_line(coocol_box2_lineA[0],coocol_box2_lineA[1],box2line_active[0],box2line_active[1],fill="red", dash = (3,3)) )
                    if coocol_box2_lineB[3]==0:
                        viz.append( canvas_1.create_line(coocol_box2_lineB[0],coocol_box2_lineB[1],box2line_active[2],box2line_active[3],fill="blue") )
                    else:
                        viz.append( canvas_1.create_line(coocol_box2_lineB[0],coocol_box2_lineB[1],box2line_active[2],box2line_active[3],fill="red", dash = (3,3)) )

                #if ((box1line_active[1]/box1line_active[3])>0.95) and ((box1line_active[1]/box1line_active[3])<1.05) or ((box1line_active[3]/box1line_active[1])>0.95) and ((box1line_active[3]/box1line_active[1])<1.05):
                #        print(1)
                #if ((box2line_active[1]/box2line_active[3])>0.95) and ((box2line_active[1]/box2line_active[3])<1.05) or ((box2line_active[3]/box2line_active[1])>0.95) and ((box2line_active[3]/box2line_active[1])<1.05):
                #        print(2)
                
                #if box1line_active[1] == box1line_active[3]:
                #        print(1)
                #if box2line_active[1] == box2line_active[3]:
                #        print(2)

                #print(coocol_box1_lineA[4])
                #print(coocol_box1_lineB[4])
                #print(coocol_box2_lineA[4])
                #print(coocol_box2_lineB[4])

                #if distpen_l1A <= r_pinball and (coocol_l1A[3]==0 or lines_edge_to_edge==1):
                    #unitvec_n = divide( [ l1.coovertex[0] - coocol_l1A[0] , l1.coovertex[1] - coocol_l1A[1] , 0 ] , distpen_l1A ) # pekar alltid mot linje 1

                if distpen_box1_lineA >0 and coocol_box1_lineA[3]==0:
                    unitvec_n = divide( [ coocol_box1_lineA[0] - box1line_active[0] , coocol_box1_lineA[1] - box1line_active[1] , 0 ] , distpen_box1_lineA )
                    #unitvec_n = coocol_box1_lineA[4]
                    #print(unitvec_n)
                    unitvec_t = Rotzvec3_90degcw(unitvec_n)
                    #print(sqrt( (box1line_active[0] - coocol_box1_lineA[0])**2 + (box1line_active[1] - coocol_box1_lineA[1])**2 ))
                    
                    ### TESTING
                    #if ((box1line_active[1]/box1line_active[3])>0.95) and ((box1line_active[1]/box1line_active[3])<1.05) or ((box1line_active[3]/box1line_active[1])>0.95) and ((box1line_active[3]/box1line_active[1])<1.05):
                    #    print(1)
                    #if ((box2line_active[1]/box2line_active[3])>0.95) and ((box2line_active[1]/box2line_active[3])<1.05) or ((box2line_active[3]/box2line_active[1])>0.95) and ((box2line_active[3]/box2line_active[1])<1.05):
                    #    print(2)
                    
                    #acoocol_box1_lineA = ClosestPointOnLineSegmentPerpDist(box1line_active[0],box1line_active[1],box2line_active,box2.coo_geom_center[0],box2.coo_geom_center[1])
                    #acoocol_box1_lineB = ClosestPointOnLineSegmentPerpDist(box1line_active[2],box1line_active[3],box2line_active,box2.coo_geom_center[0],box2.coo_geom_center[1])
                    #acoocol_box2_lineA = ClosestPointOnLineSegmentPerpDist(box2line_active[0],box2line_active[1],box1line_active,box1.coo_geom_center[0],box1.coo_geom_center[1])
                    #acoocol_box2_lineB = ClosestPointOnLineSegmentPerpDist(box2line_active[2],box2line_active[3],box1line_active,box1.coo_geom_center[0],box1.coo_geom_center[1])
                    
                    #print("distpen_box1_lineA "+str(distpen_box1_lineA))
                    #print("unitvec_n: "+str(unitvec_n))
                    #print(box2line_active[2]-box2line_active[0],box2line_active[3]-box2line_active[1])
                    #print(box2line_active)
                    ###
                    
                    #print(unitvec_n)
                    #canvas_1.create_line(box2line_active[0],box2line_active[1],box2line_active[2],box2line_active[3],fill="red")

                    Contact_line_line(box1,box2,[coocol_box1_lineA[0],coocol_box1_lineA[1],0],unitvec_n,unitvec_t)
                    
                    pendist = distpen_box1_lineA
                    #print("box1 coo1: "+str(box1.coo1))
                    box1.coo1 = add( box1.coo1 , multiply(box1.m/(box1.m+box2.m)*pendist,[unitvec_n[0],unitvec_n[1]]) )
                    #print("box1 coo1: "+str(box1.coo1))
                    box1.coo0 = add( box1.coo0 , multiply(box1.m/(box1.m+box2.m)*pendist,[unitvec_n[0],unitvec_n[1]]) )
                    box1.coom1 = add( box1.coom1 , multiply(box1.m/(box1.m+box2.m)*pendist,[unitvec_n[0],unitvec_n[1]]) )
                    box1.coovertex = add( box1.coovertex , multiply(box1.m/(box1.m+box2.m)*pendist,[unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]]) )

                    #print("box2 coo1: "+str(box2.coo1))
                    box2.coo1 = subtract( box2.coo1 , multiply(box2.m/(box2.m+box2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
                    #print("box2 coo1: "+str(box2.coo1))
                    box2.coo0 = subtract( box2.coo0 , multiply(box2.m/(box2.m+box2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
                    box2.coom1 = subtract( box2.coom1 , multiply(box2.m/(box2.m+box2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
                    box2.coovertex = subtract( box2.coovertex , multiply(box2.m/(box2.m+box2.m)*pendist, [unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]]))
                    
                '''
                ##################################### NEEDED? ###################################
                if distpen_box1_lineB >0 and coocol_box1_lineB[3]==0:
                    unitvec_n = divide( [ coocol_box1_lineB[0] - box1line_active[2] , coocol_box1_lineB[1] - box1line_active[3] , 0 ] , distpen_box1_lineB )
                    #unitvec_n = coocol_box1_lineB[4]
                    print(unitvec_n)
                    unitvec_t = Rotzvec3_90degcw(unitvec_n)
                    
                    ### TESTING
                    #if ((box1line_active[1]/box1line_active[3])>0.95) and ((box1line_active[1]/box1line_active[3])<1.05) or ((box1line_active[3]/box1line_active[1])>0.95) and ((box1line_active[3]/box1line_active[1])<1.05):
                    #    print(1)
                    #if ((box2line_active[1]/box2line_active[3])>0.95) and ((box2line_active[1]/box2line_active[3])<1.05) or ((box2line_active[3]/box2line_active[1])>0.95) and ((box2line_active[3]/box2line_active[1])<1.05):
                    #    print(2)
                    
                    #acoocol_box1_lineA = ClosestPointOnLineSegmentPerpDist(box1line_active[0],box1line_active[1],box2line_active,box2.coo_geom_center[0],box2.coo_geom_center[1])
                    #acoocol_box1_lineB = ClosestPointOnLineSegmentPerpDist(box1line_active[2],box1line_active[3],box2line_active,box2.coo_geom_center[0],box2.coo_geom_center[1])
                    #acoocol_box2_lineA = ClosestPointOnLineSegmentPerpDist(box2line_active[0],box2line_active[1],box1line_active,box1.coo_geom_center[0],box1.coo_geom_center[1])
                    #acoocol_box2_lineB = ClosestPointOnLineSegmentPerpDist(box2line_active[2],box2line_active[3],box1line_active,box1.coo_geom_center[0],box1.coo_geom_center[1])
                    
                    #print("distpen_box1_lineB "+str(distpen_box1_lineB))
                    #print("unitvec_n: "+str(unitvec_n))
                    #print(box2line_active[2]-box2line_active[0],box2line_active[3]-box2line_active[1])
                    #print(box2line_active)
                    ###

                    #print(unitvec_n)
                    #canvas_1.create_line(box2line_active[0],box2line_active[1],box2line_active[2],box2line_active[3],fill="orange",width=4)
                    #canvas_1.create_line(100,300,400,300,fill="red")

                    Contact_line_line(box1,box2,[coocol_box1_lineB[0],coocol_box1_lineB[1],0],unitvec_n,unitvec_t)
                    
                    
                    pendist = distpen_box1_lineB
                    #print("box1 coo1: "+str(box1.coo1))
                    box1.coo1 = add( box1.coo1 , multiply(box1.m/(box1.m+box2.m)*pendist,[unitvec_n[0],unitvec_n[1]]) )
                    #print("box1 coo1: "+str(box1.coo1))
                    box1.coo0 = add( box1.coo0 , multiply(box1.m/(box1.m+box2.m)*pendist,[unitvec_n[0],unitvec_n[1]]) )
                    box1.coom1 = add( box1.coom1 , multiply(box1.m/(box1.m+box2.m)*pendist,[unitvec_n[0],unitvec_n[1]]) )
                    box1.coovertex = add( box1.coovertex , multiply(box1.m/(box1.m+box2.m)*pendist,[unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]]) )

                    #print("box2 coo1: "+str(box2.coo1))
                    box2.coo1 = subtract( box2.coo1 , multiply(box2.m/(box2.m+box2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
                    #print("box2 coo1: "+str(box2.coo1))
                    box2.coo0 = subtract( box2.coo0 , multiply(box2.m/(box2.m+box2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
                    box2.coom1 = subtract( box2.coom1 , multiply(box2.m/(box2.m+box2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
                    box2.coovertex = subtract( box2.coovertex , multiply(box2.m/(box2.m+box2.m)*pendist, [unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]]))
                '''

                #################################################################################
                    
                if distpen_box2_lineA >0 and coocol_box2_lineA[3]==0:
                    unitvec_n = divide( [ box2line_active[0] - coocol_box2_lineA[0] , box2line_active[1] - coocol_box2_lineA[1] , 0 ] , distpen_box2_lineA )
                    #unitvec_n = coocol_box2_lineA[4]
                    #print(unitvec_n)
                    unitvec_t = Rotzvec3_90degcw(unitvec_n)
                    #print(sqrt( (coocol_box2_lineA[0] - box2line_active[0])**2 + (coocol_box2_lineA[1] - box2line_active[1])**2 ))
                    
                    ### TESTING
                    #if box1line_active[1] == box1line_active[3]:
                    #    print(1)
                    #if box2line_active[1] == box2line_active[3]:
                    #    print(2)
                    
                    #acoocol_box1_lineA = ClosestPointOnLineSegmentPerpDist(box1line_active[0],box1line_active[1],box2line_active,box2.coo_geom_center[0],box2.coo_geom_center[1])
                    #acoocol_box1_lineB = ClosestPointOnLineSegmentPerpDist(box1line_active[2],box1line_active[3],box2line_active,box2.coo_geom_center[0],box2.coo_geom_center[1])
                    #acoocol_box2_lineA = ClosestPointOnLineSegmentPerpDist(box2line_active[0],box2line_active[1],box1line_active,box1.coo_geom_center[0],box1.coo_geom_center[1])
                    #acoocol_box2_lineB = ClosestPointOnLineSegmentPerpDist(box2line_active[2],box2line_active[3],box1line_active,box1.coo_geom_center[0],box1.coo_geom_center[1])
                    
                    #coocol_box2_lineA = ClosestPointOnLineSegmentPerpDist(box2line_active[0],box2line_active[1],box1line_active,box1.coo_geom_center[0],box1.coo_geom_center[1])
                    #print("distpen_box2_lineA "+str(distpen_box2_lineA))
                    #print("unitvec_n: "+str(unitvec_n))
                    #print(box1line_active[2]-box1line_active[0],box1line_active[3]-box1line_active[1])
                    #print(unitvec_n)
                    ###

                    Contact_line_line(box1,box2,coocol_box2_lineA,unitvec_n,unitvec_t)
                    
                    
                    pendist = distpen_box2_lineA
                    #print("box1 coo1: "+str(box1.coo1))
                    box1.coo1 = add( box1.coo1 , multiply(box1.m/(box1.m+box2.m)*pendist,[unitvec_n[0],unitvec_n[1]]) )
                    #print("box1 coo1: "+str(box1.coo1))
                    box1.coo0 = add( box1.coo0 , multiply(box1.m/(box1.m+box2.m)*pendist,[unitvec_n[0],unitvec_n[1]]) )
                    box1.coom1 = add( box1.coom1 , multiply(box1.m/(box1.m+box2.m)*pendist,[unitvec_n[0],unitvec_n[1]]) )
                    box1.coovertex = add( box1.coovertex , multiply(box1.m/(box1.m+box2.m)*pendist,[unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]]) )

                    #print("box2 coo1: "+str(box2.coo1))
                    box2.coo1 = subtract( box2.coo1 , multiply(box2.m/(box1.m+box2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
                    #print("box2 coo1: "+str(box2.coo1))
                    box2.coo0 = subtract( box2.coo0 , multiply(box2.m/(box1.m+box2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
                    box2.coom1 = subtract( box2.coom1 , multiply(box2.m/(box1.m+box2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
                    box2.coovertex = subtract( box2.coovertex , multiply(box2.m/(box1.m+box2.m)*pendist, [unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]]))
                    
                '''
                if distpen_box2_lineB >0 and coocol_box2_lineB[3]==0:
                    unitvec_n = divide( [ box2line_active[2] - coocol_box2_lineB[0] , box2line_active[3] - coocol_box2_lineB[1] , 0 ] , distpen_box2_lineB )
                    #unitvec_n = coocol_box2_lineB[4]
                    print(unitvec_n)
                    unitvec_t = Rotzvec3_90degcw(unitvec_n)
                    #print(sqrt( (coocol_box2_lineB[0] - box2line_active[0])**2 + (coocol_box2_lineB[1] - box2line_active[1])**2 ))
                    
                    ### TESTING
                    #if box1line_active[1] == box1line_active[3]:
                    #    print(1)
                    #if box2line_active[1] == box2line_active[3]:
                    #    print(2)
                    
                    #acoocol_box1_lineA = ClosestPointOnLineSegmentPerpDist(box1line_active[0],box1line_active[1],box2line_active,box2.coo_geom_center[0],box2.coo_geom_center[1])
                    #acoocol_box1_lineB = ClosestPointOnLineSegmentPerpDist(box1line_active[2],box1line_active[3],box2line_active,box2.coo_geom_center[0],box2.coo_geom_center[1])
                    #acoocol_box2_lineA = ClosestPointOnLineSegmentPerpDist(box2line_active[0],box2line_active[1],box1line_active,box1.coo_geom_center[0],box1.coo_geom_center[1])
                    #acoocol_box2_lineB = ClosestPointOnLineSegmentPerpDist(box2line_active[2],box2line_active[3],box1line_active,box1.coo_geom_center[0],box1.coo_geom_center[1])
                    
                    #coocol_box2_lineB = ClosestPointOnLineSegmentPerpDist(box2line_active[2],box2line_active[3],box1line_active,box1.coo_geom_center[0],box1.coo_geom_center[1])
                    #print("distpen_box2_lineB "+str(distpen_box2_lineB))
                    #print("unitvec_n: "+str(unitvec_n))
                    #print(box1line_active[2]-box1line_active[0],box1line_active[3]-box1line_active[1])
                    #print(unitvec_n)
                    ###

                    Contact_line_line(box1,box2,coocol_box2_lineB,unitvec_n,unitvec_t)
                    
                    
                    pendist = distpen_box2_lineB
                    #print("box1 coo1: "+str(box1.coo1))
                    box1.coo1 = add( box1.coo1 , multiply(box1.m/(box1.m+box2.m)*pendist,[unitvec_n[0],unitvec_n[1]]) )
                    #print("box1 coo1: "+str(box1.coo1))
                    box1.coo0 = add( box1.coo0 , multiply(box1.m/(box1.m+box2.m)*pendist,[unitvec_n[0],unitvec_n[1]]) )
                    box1.coom1 = add( box1.coom1 , multiply(box1.m/(box1.m+box2.m)*pendist,[unitvec_n[0],unitvec_n[1]]) )
                    box1.coovertex = add( box1.coovertex , multiply(box1.m/(box1.m+box2.m)*pendist,[unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]]) )

                    #print("box2 coo1: "+str(box2.coo1))
                    box2.coo1 = subtract( box2.coo1 , multiply(box2.m/(box1.m+box2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
                    #print("box2 coo1: "+str(box2.coo1))
                    box2.coo0 = subtract( box2.coo0 , multiply(box2.m/(box1.m+box2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
                    box2.coom1 = subtract( box2.coom1 , multiply(box2.m/(box1.m+box2.m)*pendist, [unitvec_n[0],unitvec_n[1]] ))
                    box2.coovertex = subtract( box2.coovertex , multiply(box2.m/(box1.m+box2.m)*pendist, [unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]] ))
                '''

def Collision_Polygon_FixedLine(box,fl):
    
    r_pinball = fl.hhalf

    #global colcounter

    if sqrt( (box.coo0[0]-fl.coo0[0])**2 + (box.coo0[1]-fl.coo0[1])**2 ) <= (box.r + fl.r + r_pinball):

        #indices = [ [0,1,2,3] , [2,3,4,5] , [4,5,6,7] , [6,7,0,1] ]

        #for i in range(0,len(indices)):
        for i in range(0,box.n_linesegments):
            #boxline_active = [box.coovertex[ indices[i][0] ],box.coovertex[ indices[i][1] ],box.coovertex[ indices[i][2] ],box.coovertex[ indices[i][3] ]]
            if i == (box.n_linesegments-1): # if i is at the last xy-pair in the list
                boxline_active = [box.coovertex[2*i],box.coovertex[2*i+1],box.coovertex[0],box.coovertex[1]]
            else: # else, count as usual
                boxline_active = [box.coovertex[2*i],box.coovertex[2*i+1],box.coovertex[2*i+2],box.coovertex[2*i+3]]


            coocol_flA = ClosestPointOnLineSegmentAbsDist(fl.coovertex[0],fl.coovertex[1],boxline_active)
            coocol_flB = ClosestPointOnLineSegmentAbsDist(fl.coovertex[2],fl.coovertex[3],boxline_active)
            #distint_flA = sqrt( (coocol_flA[0]-fl.coovertex[0])**2 + (coocol_flA[1]-fl.coovertex[1])**2 )
            #distint_flB = sqrt( (coocol_flB[0]-fl.coovertex[2])**2 + (coocol_flB[1]-fl.coovertex[3])**2 )

            distint_flA = coocol_flA[2]
            distint_flB = coocol_flB[2]

            coocol_boxA = ClosestPointOnLineSegmentAbsDist(boxline_active[0],boxline_active[1],fl.coovertex)
            coocol_boxB = ClosestPointOnLineSegmentAbsDist(boxline_active[2],boxline_active[3],fl.coovertex)
            #distint_boxA = sqrt( (coocol_boxA[0]-boxline_active[0])**2 + (coocol_boxA[1]-boxline_active[1])**2 )
            #distint_boxB = sqrt( (coocol_boxB[0]-boxline_active[2])**2 + (coocol_boxB[1]-boxline_active[3])**2 )
            
            distint_boxA = coocol_boxA[2]
            distint_boxB = coocol_boxB[2]

            if debug_col_lines==1:
                viz.append( canvas_1.create_line(coocol_boxA[0],coocol_boxA[1],boxline_active[0],boxline_active[1],fill="green") )
                viz.append( canvas_1.create_line(coocol_boxB[0],coocol_boxB[1],boxline_active[2],boxline_active[3],fill="green") )
                viz.append( canvas_1.create_line(coocol_flA[0],coocol_flA[1],fl.coovertex[0],fl.coovertex[1],fill="green") )
                viz.append( canvas_1.create_line(coocol_flB[0],coocol_flB[1],fl.coovertex[2],fl.coovertex[3],fill="green") )

            if distint_boxA <= r_pinball and coocol_boxA[3]==0:
                unitvec_n = divide( [ boxline_active[0] - coocol_boxA[0] , boxline_active[1] - coocol_boxA[1] , 0 ] , distint_boxA )
                unitvec_t = Rotzvec3_90degcw(unitvec_n)

                coocol_hcomp = add( [coocol_boxA[0],coocol_boxA[1],0], multiply( fl.hhalf, unitvec_n ) )

                Contact_dyn_statline(box,fl,coocol_hcomp,unitvec_n,unitvec_t)
        
                box.coo1 = box.coo1 + multiply([unitvec_n[0],unitvec_n[1]],(r_pinball-distint_boxA))
                box.coo0 = box.coo0 + multiply([unitvec_n[0],unitvec_n[1]],(r_pinball-distint_boxA))
                box.coom1 = box.coom1 + multiply([unitvec_n[0],unitvec_n[1]],(r_pinball-distint_boxA))
                box.coom_geom_center = box.coo_geom_center + multiply([unitvec_n[0],unitvec_n[1]],(r_pinball-distint_boxA))

                #box.coovertex = box.coovertex + multiply([unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]],(r_pinball-distint_boxA))
                for j in range(box.n_linesegments):
                    box.coovertex[2*j] = box.coovertex[2*j] + unitvec_n[0]*(r_pinball-distint_boxA)
                    box.coovertex[2*j+1] = box.coovertex[2*j+1] + unitvec_n[1]*(r_pinball-distint_boxA)

                #colcounter = colcounter + 1
                #print(str(colcounter)+" - boxline "+str(i)+", edge A, indicator "+str(coocol_boxA[3]))
            '''
            if distint_boxB <= r_pinball and coocol_boxB[2]==0:
                unitvec_n = divide( [ boxline_active[2] - coocol_boxB[0] , boxline_active[3] - coocol_boxB[1] , 0 ] , distint_boxB )
                unitvec_t = Rotzvec3_90degcw(unitvec_n)

                coocol_hcomp = add( [coocol_boxB[0],coocol_boxB[1],0], multiply( fl.hhalf, unitvec_n ) )

                Contact_dyn_statline(box,fl,coocol_hcomp,unitvec_n,unitvec_t)
            
                box.coo1 = box.coo1 + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_boxB))
                box.coo0 = box.coo0 + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_boxB))
                box.coom1 = box.coom1 + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_boxB))
                #box.cooA = box.cooA + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_boxB))
                #box.cooB = box.cooB + multiply([unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_boxB))
                box.coovertex = box.coovertex + multiply([unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]],fac*(r_pinball-distint_boxB))

                #colcounter = colcounter + 1
                #print(str(colcounter)+" - boxline "+str(i)+", edge B, indicator "+str(coocol_boxB[3]))
            '''
            if distint_flA <= r_pinball and coocol_flA[3]==0:
                unitvec_n = divide( [ coocol_flA[0] - fl.coovertex[0] , coocol_flA[1] - fl.coovertex[1] , 0 ] , distint_flA )
                unitvec_t = Rotzvec3_90degcw(unitvec_n)

                Contact_dyn_statline(box,fl,coocol_flA,unitvec_n,unitvec_t)

                box.coo1 = box.coo1 + multiply([unitvec_n[0],unitvec_n[1]],(r_pinball-distint_flA))
                box.coo0 = box.coo0 + multiply([unitvec_n[0],unitvec_n[1]],(r_pinball-distint_flA))
                box.coom1 = box.coom1 + multiply([unitvec_n[0],unitvec_n[1]],(r_pinball-distint_flA))
                box.coom_geom_center = box.coo_geom_center + multiply([unitvec_n[0],unitvec_n[1]],(r_pinball-distint_flA))
                
                #box.coovertex = box.coovertex + multiply([unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]],(r_pinball-distint_flA))
                for j in range(box.n_linesegments):
                    box.coovertex[2*j] = box.coovertex[2*j] + unitvec_n[0]*(r_pinball-distint_flA)
                    box.coovertex[2*j+1] = box.coovertex[2*j+1] + unitvec_n[1]*(r_pinball-distint_flA)

                #colcounter = colcounter + 1
                #print(str(colcounter)+" - flA - "+str(coocol_flA[3]))
        
            if distint_flB <= r_pinball and coocol_flB[3]==0:
                unitvec_n = divide( [ coocol_flB[0] - fl.coovertex[2] , coocol_flB[1] - fl.coovertex[3] , 0 ] , distint_flB )
                unitvec_t = Rotzvec3_90degcw(unitvec_n)

                Contact_dyn_statline(box,fl,coocol_flB,unitvec_n,unitvec_t)
                
                box.coo1 = box.coo1 + multiply([unitvec_n[0],unitvec_n[1]],(r_pinball-distint_flB))
                box.coo0 = box.coo0 + multiply([unitvec_n[0],unitvec_n[1]],(r_pinball-distint_flB))
                box.coom1 = box.coom1 + multiply([unitvec_n[0],unitvec_n[1]],(r_pinball-distint_flB))
                box.coom_geom_center = box.coo_geom_center + multiply([unitvec_n[0],unitvec_n[1]],(r_pinball-distint_flB))
                
                #box.coovertex = box.coovertex + multiply([unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1],unitvec_n[0],unitvec_n[1]],(r_pinball-distint_flB))
                for j in range(box.n_linesegments):
                    box.coovertex[2*j] = box.coovertex[2*j] + unitvec_n[0]*(r_pinball-distint_flB)
                    box.coovertex[2*j+1] = box.coovertex[2*j+1] + unitvec_n[1]*(r_pinball-distint_flB)

                #colcounter = colcounter + 1
                #print(str(colcounter)+" - flB - "+str(coocol_flB[3]))
                
            #if obj[10] is fl:
            #    #print(coocol_boxA)
            #    Dp.moveto( coocol_boxA[0],coocol_boxA[1] )

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
        self.F_D_form = array([0.0, 0.0])
        self.F_R = array([0.0, 0.0])
        # Draw
        self.canvas = canvas
        self.ball = canvas.create_oval(x-self.r,y-self.r,x+self.r,y+self.r,fill="white")
        self.thetaline = canvas.create_line(self.coo1[0],self.coo1[1],self.coothetaline[0],self.coothetaline[1])
        #self.arrow_F = canvas.create_line(self.coo1[0],self.coo1[1],self.F_other[0],self.F_other[1],arrow=LAST,fill="grey")
        #self.text_F = canvas_1.create_text(self.coo0[0],self.coo0[1],text=str(self.F_other_abs),font=("arial",8),fill="grey")
        
        self.arrow_g = [ sf_farrow1*tanh((self.m*grav[0])/sf_farrow2) , sf_farrow1*tanh((self.m*grav[1])/sf_farrow2) ]
        self.arrow_F_D_form = [ (sf_farrow1*tanh(self.F_D_form[0]/sf_farrow2)) , (sf_farrow1*tanh(self.F_D_form[1]/sf_farrow2)) ]
        self.arrow_g_draw = canvas.create_line(self.coo1[0],self.coo1[1],self.coo1[0]+self.arrow_g[0],self.coo1[1]+self.arrow_g[1],arrow=LAST,fill="red")
        #self.offset_text_coo = -20
        #self.text_coo = canvas_1.create_text(self.coo0[0],self.coo0[1]+self.offset_text_coo,text="("+str(round(self.coo0[0]))+", "+str(round(self.coo0[1]))+")",font=("arial",8))        
        self.arrow_F_D_draw = canvas.create_line(self.coo1[0],self.coo1[1],self.coo1[0]-self.arrow_F_D_form[0],self.coo1[1]-self.arrow_F_D_form[1],arrow=LAST,fill="blue")
        

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
        self.arrow_F_D_form = [ (sf_farrow1*tanh(self.F_D_form[0]/sf_farrow2)) , (sf_farrow1*tanh(self.F_D_form[1]/sf_farrow2)) ]
        self.canvas.coords(self.arrow_g_draw,self.coo1[0],self.coo1[1],self.coo1[0]+self.arrow_g[0],self.coo1[1]+self.arrow_g[1])
        self.canvas.coords(self.arrow_F_D_draw,self.coo1[0],self.coo1[1],self.coo1[0]-self.arrow_F_D_form[0],self.coo1[1]-self.arrow_F_D_form[1])
        

        #self.canvas.coords(self.text_coo,self.coo0[0],self.coo0[1]+self.offset_text_coo)
        #self.canvas.itemconfigure(self.text_coo,text="("+str(round(self.coo0[0]))+", "+str(round(self.coo0[1]))+")")

    def integrate(self):

        self.vel = subtract(self.coo1,self.coo0)/dt
        self.absvel = sqrt( (self.vel[0])**2 + (self.vel[1])**2 )
        self.unitvec_vel = [Safediv(self.vel[0],self.absvel),Safediv(self.vel[1],self.absvel)]
        self.A_ABfacingdirection = pi*self.r**2

        # Drag force in the chosen medium (air, water, etc.)
        self.F_D_form = multiply( 0.5*rho_medium*0.50*self.A_ABfacingdirection , multiply( [ (self.vel[0])**2 , (self.vel[1])**2 ] , self.unitvec_vel ) )
        F_crit = self.m*self.absvel/(2*dt)
        if sqrt( (self.F_D_form[0])**2 + (self.F_D_form[1])**2 ) >= F_crit:
            self.F_D_form = multiply(F_crit,self.unitvec_vel)
            #print("                                      Drag force limit exceeded" + str(self.F_D))
        self.F_other = subtract( self.F_other , self.F_D_form )
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

class Object_Polygon:
    def __init__(self,canvas,x,y,pointsvector,t,rho,C_el,C_fric,ecc_u,ecc_v):
        # General
        self.canvas = canvas
        # Geometrical properties
        #self.w = w
        #self.h = h
        self.t = t
        self.C_fric = C_fric
        self.C_el = C_el
        # Coordinates
        self.coo1 = [ x, y ]
        #self.k = array([ Safediv((self.coo1[5]-self.coo1[3]),(self.coo1[4]-self.coo1[2])) , Safediv((self.coo1[7]-self.coo1[5]),(self.coo1[6]-self.coo1[4])) , Safediv((self.coo1[9]-self.coo1[7]),(self.coo1[8]-self.coo1[6])) , Safediv((self.coo1[3]-self.coo1[9]),(self.coo1[2]-self.coo1[8])) ])
        #print(self.k)
        #self.coo1 = [x,y]
        self.coo0 = self.coo1
        self.coom1 = self.coo1
        self.theta1 = 0
        self.theta0 = 0
        self.thetam1 = self.theta1
        self.ecc_u = ecc_u
        self.ecc_v = ecc_v
        self.n_linesegments = int(len(pointsvector)*0.5)

        #self.coovertexloc = [ -0.5*self.w+self.ecc_u*self.w , -0.5*self.h+self.ecc_v*self.h , 0.5*self.w+self.ecc_u*self.w , -0.5*self.h+self.ecc_v*self.h , 0.5*self.w+self.ecc_u*self.w , 0.5*self.h+self.ecc_v*self.h , -0.5*self.w+self.ecc_u*self.w , 0.5*self.h+self.ecc_v*self.h ]
        self.coovertexloc = [None]*2*self.n_linesegments
        for i in range(self.n_linesegments):
            self.coovertexloc[2*i] = pointsvector[2*i]+self.ecc_u
            self.coovertexloc[2*i+1] = pointsvector[2*i+1]+self.ecc_v

        #self.coovertexloc = [ -0.5*w , -0.5*h , 0.5*w , -0.5*h , 0.5*w , 0.5*h , -0.5*w , 0.5*h ]
        self.coovertexlocrot = Rotxypairsinvec(self.coovertexloc,self.theta1)
        #self.coovertex = [ self.coo1[0]+self.coovertexlocrot[0] , self.coo1[1]+self.coovertexlocrot[1] , self.coo1[0]+self.coovertexlocrot[2] , self.coo1[1]+self.coovertexlocrot[3] , self.coo1[0]+self.coovertexlocrot[4] , self.coo1[1]+self.coovertexlocrot[5] , self.coo1[0]+self.coovertexlocrot[6] , self.coo1[1]+self.coovertexlocrot[7] ]
        
        self.coovertex = [None]*2*self.n_linesegments
        self.r = 0
        self.coo_geom_center = [0,0]
        for i in range(self.n_linesegments):
            # calculate vertex coordinates
            self.coovertex[2*i] = self.coo1[0]+self.coovertexlocrot[2*i]
            self.coovertex[2*i+1] = self.coo1[1]+self.coovertexlocrot[2*i+1]
            # calculate point furthest away from center of mass (for contact search)
            rnew = sqrt( (self.coovertexloc[2*i])**2 + (self.coovertexloc[2*i+1])**2 )
            if rnew > self.r:
                self.r = rnew

            # also geometrical center of polygon while looping
            self.coo_geom_center[0] = self.coo_geom_center[0] + self.coovertex[2*i]
            self.coo_geom_center[1] = self.coo_geom_center[1] + self.coovertex[2*i+1]
        self.coo_geom_center = divide(self.coo_geom_center,self.n_linesegments)
            




        ## calculate point furthest away from center of mass (for contact search)
        #com_to_A = sqrt( (self.coovertexloc[0])**2 + (self.coovertexloc[1])**2 )
        #com_to_B = sqrt( (self.coovertexloc[2])**2 + (self.coovertexloc[3])**2 )
        #com_to_C = sqrt( (self.coovertexloc[4])**2 + (self.coovertexloc[5])**2 )
        #com_to_D = sqrt( (self.coovertexloc[6])**2 + (self.coovertexloc[7])**2 )
        #self.r = max(com_to_A,com_to_B,com_to_C,com_to_D)

        # geometrical center of polygon
        #self.coo_geom_center = ( (self.coovertex[0]+self.coovertex[2]+self.coovertex[4]+self.coovertex[6])/4 , (self.coovertex[1]+self.coovertex[3]+self.coovertex[5]+self.coovertex[7])/4 )
        
        
        # calculate boundary normals and angles at no rotation
        self.normals_local = [None]*2*self.n_linesegments
        self.tangents_local = [None]*2*self.n_linesegments
        self.angles_local = [None]*self.n_linesegments
        
        # calculate local normals, local tangents, and local angles
        for i in range(self.n_linesegments):
            if i == (self.n_linesegments-1): # if i is at the last xy-pair in the list
                #line = [self.coovertex[2*i],self.coovertex[2*i+1],self.coovertex[0],self.coovertex[1]]
                unitvec_t = [ self.coovertexloc[0] - self.coovertexloc[2*i] , self.coovertexloc[1] - self.coovertexloc[2*i+1] ] / sqrt( (self.coovertexloc[0] - self.coovertexloc[2*i])**2 + (self.coovertexloc[1] - self.coovertexloc[2*i+1])**2 )
                unitvec_n = Rotvec2_90degcountercw(unitvec_t)
                self.normals_local[2*i] = unitvec_n[0]
                self.normals_local[2*i+1] = unitvec_n[1]
                self.tangents_local[2*i] = unitvec_t[0]
                self.tangents_local[2*i+1] = unitvec_t[1]
                self.angles_local[i] = arctan(Safediv( (self.coovertexloc[1] - self.coovertexloc[2*i+1]) , (self.coovertexloc[0] - self.coovertexloc[2*i]) ))
            else: # else, count as usual
                #line = [self.coovertex[2*i],self.coovertex[2*i+1],self.coovertex[2*i+2],self.coovertex[2*i+3]]
                unitvec_t = [ self.coovertexloc[2*i+2] - self.coovertexloc[2*i] , self.coovertexloc[2*i+3] - self.coovertexloc[2*i+1] ] / sqrt( (self.coovertexloc[2*i+2] - self.coovertexloc[2*i])**2 + (self.coovertexloc[2*i+3] - self.coovertexloc[2*i+1])**2 )
                unitvec_n = Rotvec2_90degcountercw(unitvec_t)
                self.normals_local[2*i] = unitvec_n[0]
                self.normals_local[2*i+1] = unitvec_n[1] 
                self.tangents_local[2*i] = unitvec_t[0]
                self.tangents_local[2*i+1] = unitvec_t[1]
                self.angles_local[i] = arctan(Safediv( (self.coovertexloc[2*i+3] - self.coovertexloc[2*i+1]) , (self.coovertexloc[2*i+2] - self.coovertexloc[2*i]) ))
            #canvas_1.create_line( self.coo1[0]+unitvec_n[0]*25,self.coo1[1]+unitvec_n[1]*25,self.coo1[0]+unitvec_n[0]*50,self.coo1[1]+unitvec_n[1]*50,arrow=LAST )
        # rotate them to the starting rotation
        self.normals_local_rotated = Rotxypairsinvec(self.normals_local,self.theta1)
        self.tangents_local_rotated = Rotxypairsinvec(self.tangents_local,self.theta1)
        self.angles_local_rotated = add(self.angles_local,self.theta1)
        
        # Forces and torques
        self.F_g = [0.0,0.0]
        self.F_other = [0.0,0.0]
        self.tau_other = 0.0
        # Inertia
        #self.m = w*h*t*rho
        self.m = 40000*rho                    #temporary
        self.invm = 1/self.m
        #self.I = 0.08333*self.m*(w*w+h*h)
        self.I = 0.08333*self.m*(4100)        #temporary
        self.invImat = [ [1/self.I,0,0] , [0,0,0] , [0,0,1/self.I] ]
        # Draw
        #self.lineAB = self.canvas.create_line( self.coovertex[0],self.coovertex[1],self.coovertex[2],self.coovertex[3] )
        #self.lineBC = self.canvas.create_line( self.coovertex[2],self.coovertex[3],self.coovertex[4],self.coovertex[5] )
        #self.lineCD = self.canvas.create_line( self.coovertex[4],self.coovertex[5],self.coovertex[6],self.coovertex[7] )
        #self.lineDA = self.canvas.create_line( self.coovertex[6],self.coovertex[7],self.coovertex[0],self.coovertex[1] )
        self.lines = [None]*self.n_linesegments
        for i in range(self.n_linesegments):
            if i == (self.n_linesegments-1): # if i is at the last xy-pair in the list
                self.lines[i] = self.canvas.create_line( self.coovertex[2*i],self.coovertex[2*i+1],self.coovertex[0],self.coovertex[1] )
            else:
                self.lines[i] = self.canvas.create_line( self.coovertex[2*i],self.coovertex[2*i+1],self.coovertex[2*i+2],self.coovertex[2*i+3] )
        self.ball = self.canvas.create_oval( self.coo1[0]-3,self.coo1[1]-3,self.coo1[0]+3,self.coo1[1]+3,outline="grey",fill="grey" )

    def linesegment(self,n):
        if n == (self.n_linesegments-1): # if i is at the last xy-pair in the list
            n_output = [self.coovertex[2*n],self.coovertex[2*n+1],self.coovertex[0],self.coovertex[1]]
        else: # else, count as usual
            n_output = [self.coovertex[2*n],self.coovertex[2*n+1],self.coovertex[2*n+2],self.coovertex[2*n+3]]
        return n_output

    def draw(self):
        #self.canvas.coords( self.lineAB,self.coovertex[0],self.coovertex[1],self.coovertex[2],self.coovertex[3] )
        #self.canvas.coords( self.lineBC,self.coovertex[2],self.coovertex[3],self.coovertex[4],self.coovertex[5] )
        #self.canvas.coords( self.lineCD,self.coovertex[4],self.coovertex[5],self.coovertex[6],self.coovertex[7] )
        #self.canvas.coords( self.lineDA,self.coovertex[6],self.coovertex[7],self.coovertex[0],self.coovertex[1] )
        for i in range(self.n_linesegments):

            #canvas_1.create_text(10,20,text="dt = "+str(dt),font=("arial",8),anchor=SW)
            #viz.append( canvas_1.create_rectangle(ofs-100,ofs-100,ofs+100,ofs+100,fill="white") )
            viz.append( canvas_1.create_text(self.coovertex[2*i],self.coovertex[2*i+1],text=str(i),font=("arial",8),fill="red") )

            if i == (self.n_linesegments-1): # if i is at the last xy-pair in the list
                self.canvas.coords( self.lines[i],self.coovertex[2*i],self.coovertex[2*i+1],self.coovertex[0],self.coovertex[1] )
            else:
                self.canvas.coords( self.lines[i],self.coovertex[2*i],self.coovertex[2*i+1],self.coovertex[2*i+2],self.coovertex[2*i+3] )
        self.canvas.coords( self.ball,self.coo1[0]-3,self.coo1[1]-3,self.coo1[0]+3,self.coo1[1]+3 )



    def integrate(self):

        self.F_g = multiply(self.m,grav)
        self.F_other = add( self.F_other, self.F_g )

        self.coom1 = self.coo0
        self.coo0 = self.coo1
        self.thetam1 = self.theta0
        self.theta0 = self.theta1
        
        self.v0 = array([0,0])#divide(subtract(self.coo0,self.coo1), dt)
        self.coo1 = [ 2*self.coo0[0] - self.coom1[0] + (self.F_other[0]/self.m)*dt**2 , 2*self.coo0[1] - self.coom1[1] + (self.F_other[1]/self.m)*dt**2 ]
        self.theta1 = 2*self.theta0 - self.thetam1 + (self.tau_other/self.I)*dt**2
        self.coovertexlocrot = Rotxypairsinvec(self.coovertexloc,self.theta1)
        #self.coovertex = [ self.coo1[0]+self.coovertexlocrot[0] , self.coo1[1]+self.coovertexlocrot[1] , self.coo1[0]+self.coovertexlocrot[2] , self.coo1[1]+self.coovertexlocrot[3] , self.coo1[0]+self.coovertexlocrot[4] , self.coo1[1]+self.coovertexlocrot[5] , self.coo1[0]+self.coovertexlocrot[6] , self.coo1[1]+self.coovertexlocrot[7] ]
        
        #self.coo_geom_center = ( (self.coovertex[0]+self.coovertex[2]+self.coovertex[4]+self.coovertex[6])/4 , (self.coovertex[1]+self.coovertex[3]+self.coovertex[5]+self.coovertex[7])/4 )
        for i in range(self.n_linesegments):
            # calculate vertex coordinates
            self.coovertex[2*i] = self.coo1[0]+self.coovertexlocrot[2*i]
            self.coovertex[2*i+1] = self.coo1[1]+self.coovertexlocrot[2*i+1]
            # geometrical center of polygon
            self.coo_geom_center[0] = self.coo_geom_center[0] + self.coovertex[2*i]
            self.coo_geom_center[1] = self.coo_geom_center[1] + self.coovertex[2*i+1]
        self.coo_geom_center = divide(self.coo_geom_center,self.n_linesegments)

        self.normals_local_rotated = Rotxypairsinvec(self.normals_local,self.theta1)
        self.tangents_local_rotated = Rotxypairsinvec(self.tangents_local,self.theta1)
        self.angles_local_rotated = add(self.angles_local,self.theta1)

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
        self.coovertex = [ xA, yA, xB, yB ]
        self.coo0 = [(xA+xB)/2, (yA+yB)/2]
        # Bättre att beräkna dessa en gång, istället för varje gång under kontaktsökningen
        #self.coomax = [ max(self.cooA[0],self.cooB[0]) , max(self.cooA[1],self.cooB[1]) ]
        #self.coomin = [ min(self.cooA[0],self.cooB[0]) , min(self.cooA[1],self.cooB[1]) ]
        self.AB = sqrt( (xB-xA)**2 + (yB-yA)**2 )
        self.r = self.AB/2
        '''
        self.theta1 = arcsin( Safediv( abs(yB-yA) , self.AB ) )
        # due to the angle being calculated with arcsin, and due to sins periodicity, this is necessary. otherwise, angles may be wrong, depending on endpoint placement
        if xB<xA:
            self.theta1 = pi - self.theta1
        if yB<yA:
            self.theta1 = - self.theta1
        '''
        self.theta1 = arctan( Safediv( yB-yA , xB-xA ) )
        #self.text_test = canvas_1.create_text(self.coo0[0],self.coo0[1]-20,text=str(round(180/pi*self.theta1,3)),font=("arial",8))

        self.coovertexloc = [ -self.r , 0 , self.r , 0 ]
        self.coovertexlocrot = Rotxypairsinvec(self.coovertexloc,self.theta1)
        self.coovertex = [ self.coo0[0]+self.coovertexlocrot[0] , self.coo0[1]+self.coovertexlocrot[1] , self.coo0[0]+self.coovertexlocrot[2] , self.coo0[1]+self.coovertexlocrot[3] ]

        #self.unitvec_t = Safediv( [ (xB-xA) , (yB-yA) , 0] ,self.AB)
        #self.unitvec_n = Rotzvec3_90degcw(self.unitvec_t)

        #self.ang = arccos( (xB-xA)/self.AB )
        #self.coomidp = [ self.cooA[0] + 0.5*(self.cooB[0]-self.cooA[0]) , self.cooA[1] + 0.5*(self.cooB[1]-self.cooA[1]) ]
        #self.unitvec_n = Rotvec2( [1,0] , self.theta1-pi/2 )
        #self.normalvec = add( Rotvec2( [20,0] , self.ang-pi/2 ) , self.coomidp )
        #self.nvec = [0,0]
        #self.k = Safediv( (self.cooB[1]-self.cooA[1]) , (self.cooB[0]-self.cooA[0]) )
        #self.coocol = [0,0]
        #self.distint = 0
        #self.text_distint = canvas_1.create_text(0,0,text="",font=("arial",8))

        # Forces
        self.F_other = [0.0,0.0]
        # Draw
        self.canvas = canvas
        self.lineAB = canvas.create_line(self.coovertex[0],self.coovertex[1],self.coovertex[2],self.coovertex[3],width=h)
        self.ballA = canvas.create_oval(self.coovertex[0]-self.hhalf,self.coovertex[1]-self.hhalf,self.coovertex[0]+self.hhalf,self.coovertex[1]+self.hhalf,outline="black",fill="black")
        self.ballB = canvas.create_oval(self.coovertex[2]-self.hhalf,self.coovertex[3]-self.hhalf,self.coovertex[2]+self.hhalf,self.coovertex[3]+self.hhalf,outline="black",fill="black")
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
        self.theta1 = arctan( Safediv( yB-yA , xB-xA ) )
        
        self.coovertexloc = [ -self.r , 0 , self.r , 0 ]
        self.coovertexlocrot = Rotxypairsinvec(self.coovertexloc,self.theta1)
        self.coovertex = [ self.coo1[0]+self.coovertexlocrot[0] , self.coo1[1]+self.coovertexlocrot[1] , self.coo1[0]+self.coovertexlocrot[2] , self.coo1[1]+self.coovertexlocrot[3] ]
        
        #self.unitvec_t = Safediv([ (xB-xA) , (yB-yA) , 0 ],self.AB)
        #self.unitvec_n = Rotzvec3_90degcw(self.unitvec_t)

        self.vel = subtract(self.coo1,self.coo0)/dt
        self.absvel = sqrt( (self.vel[0])**2 + (self.vel[1])**2 )
        self.unitvec_vel = [Safediv(self.vel[0],self.absvel),Safediv(self.vel[1],self.absvel)]
        self.unitvec_velperp = Rotvec2_90degcw(self.unitvec_vel)
        self.vec_ABfacingdirection = Vecproj_vec2(self.vec_AB,self.unitvec_velperp)
        self.A_ABfacingdirection = self.t*sqrt( (self.vec_ABfacingdirection[0])**2 + (self.vec_ABfacingdirection[1])**2 )


        self.theta0 = self.theta1
        self.thetam1 = self.theta1
        self.F_other = [0,0]
        self.F_D_form = [0,0]
        self.F_R = array([0.0, 0.0]) # Behövs verkligen denna?
        self.tau_other = 0

        self.V = self.AB*h*t
        self.m = rho*self.V
        self.invm = 1/self.m
        self.I = 0.08333*self.m*(h*h+self.AB**2)
        self.invImat = [ [1/self.I,0,0] , [0,0,0] , [0,0,1/self.I] ]
        
        # Forces
        self.F_other = [0.0,0.0]
        # Draw
        self.canvas = canvas
        self.lineAB = canvas.create_line(self.coovertex[0],self.coovertex[1],self.coovertex[2],self.coovertex[3],width=h,fill="black")
        self.ballA = canvas.create_oval(self.coovertex[0]-self.hhalf,self.coovertex[1]-self.hhalf,self.coovertex[0]+self.hhalf,self.coovertex[1]+self.hhalf,outline="black",fill="black")
        self.ballB = canvas.create_oval(self.coovertex[2]-self.hhalf,self.coovertex[3]-self.hhalf,self.coovertex[2]+self.hhalf,self.coovertex[3]+self.hhalf,outline="black",fill="black")
        
        self.arrow_g = [ sf_farrow1*tanh((self.m*grav[0])/sf_farrow2) , sf_farrow1*tanh((self.m*grav[1])/sf_farrow2) ]
        self.arrow_F_D_form = [ (sf_farrow1*tanh(self.F_D_form[0]/sf_farrow2)) , (sf_farrow1*tanh(self.F_D_form[1]/sf_farrow2)) ]
        self.arrow_g_draw = canvas.create_line(self.coo1[0],self.coo1[1],self.coo1[0]+self.arrow_g[0],self.coo1[1]+self.arrow_g[1],arrow=LAST,fill="red")
        #self.offset_text_coo = -20
        #self.text_coo = canvas_1.create_text(self.coo0[0],self.coo0[1]+self.offset_text_coo,text="("+str(round(self.coo0[0]))+", "+str(round(self.coo0[1]))+")",font=("arial",8))        
        self.arrow_F_D_draw = canvas.create_line(self.coo1[0],self.coo1[1],self.coo1[0]-self.arrow_F_D_form[0],self.coo1[1]-self.arrow_F_D_form[1],arrow=LAST,fill="blue")
        
        #self.lineN = canvas.create_line(self.coomidp[0],self.coomidp[1],self.normalvec[0],self.normalvec[1],arrow=LAST)
        #self.linecol = canvas.create_line( 0,0,1,0, dash = (2,2) )
        



    def integrate(self):
        
        self.vec_AB = [ (self.coovertex[2]-self.coovertex[0]) , (self.coovertex[3]-self.coovertex[1]) ]
        # Ersätt med en Safediv för vektorer. Nämnaren är 0 vid programstart.
        self.vel = subtract(self.coo1,self.coo0)/dt
        self.absvel = sqrt( (self.vel[0])**2 + (self.vel[1])**2 )
        self.unitvec_vel = [Safediv(self.vel[0],self.absvel),Safediv(self.vel[1],self.absvel)]
        self.unitvec_velperp = Rotvec2_90degcw(self.unitvec_vel)
        self.vec_ABfacingdirection = Vecproj_vec2(self.vec_AB,self.unitvec_velperp)
        self.A_ABfacingdirection = self.t*sqrt( (self.vec_ABfacingdirection[0])**2 + (self.vec_ABfacingdirection[1])**2 )

        # Drag force in the chosen medium (air, water, etc.)
        self.F_D_form = multiply( 0.5*rho_medium*1.98*self.A_ABfacingdirection , multiply( [ (self.vel[0])**2 , (self.vel[1])**2 ] , self.unitvec_vel ) )
        F_crit = self.m*self.absvel/(2*dt)
        if sqrt( (self.F_D_form[0])**2 + (self.F_D_form[1])**2 ) >= F_crit:
            self.F_D_form = multiply(F_crit,self.unitvec_vel)
            #print("                                      Drag force limit exceeded" + str(self.F_D))
        self.F_other = subtract( self.F_other , self.F_D_form )
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
        #self.k = Safediv( (self.cooB[1]-self.cooA[1]) , (self.cooB[0]-self.cooA[0]) )
        #self.cooAloc = Rotvec2([-0.5*self.AB,0],self.theta1)
        #self.cooBloc = Rotvec2([0.5*self.AB,0],self.theta1)
        #self.cooA = add(self.coo1 , self.cooAloc)
        #self.cooB = add(self.coo1 , self.cooBloc)
        #self.coomax = [ max(self.cooA[0],self.cooB[0]) , max(self.cooA[1],self.cooB[1]) ] #ta bort?
        #self.coomin = [ min(self.cooA[0],self.cooB[0]) , min(self.cooA[1],self.cooB[1]) ] #ta bort?

        #print(self.cooA)
        #print(self.cooB)
        
        self.coovertexlocrot = Rotxypairsinvec(self.coovertexloc,self.theta1)
        self.coovertex = [ self.coo1[0]+self.coovertexlocrot[0] , self.coo1[1]+self.coovertexlocrot[1] , self.coo1[0]+self.coovertexlocrot[2] , self.coo1[1]+self.coovertexlocrot[3] ]
        
        #self.unitvec_t = Safediv([ (self.coovertex[2]-self.coovertex[0]) , (self.coovertex[3]-self.coovertex[1]) , 0 ],self.AB)
        #self.unitvec_n = Rotzvec3_90degcw(self.unitvec_t)

        self.canvas.coords(self.lineAB,self.coovertex[0],self.coovertex[1],self.coovertex[2],self.coovertex[3])
        self.canvas.coords(self.ballA,self.coovertex[0]-self.hhalf,self.coovertex[1]-self.hhalf,self.coovertex[0]+self.hhalf,self.coovertex[1]+self.hhalf)
        self.canvas.coords(self.ballB,self.coovertex[2]-self.hhalf,self.coovertex[3]-self.hhalf,self.coovertex[2]+self.hhalf,self.coovertex[3]+self.hhalf)

        self.arrow_g = [ sf_farrow1*tanh((self.m*grav[0])/sf_farrow2) , sf_farrow1*tanh((self.m*grav[1])/sf_farrow2) ]
        self.arrow_F_D_form = [ (sf_farrow1*tanh(self.F_D_form[0]/sf_farrow2)) , (sf_farrow1*tanh(self.F_D_form[1]/sf_farrow2)) ]
        self.canvas.coords(self.arrow_g_draw,self.coo1[0],self.coo1[1],self.coo1[0]+self.arrow_g[0],self.coo1[1]+self.arrow_g[1])
        self.canvas.coords(self.arrow_F_D_draw,self.coo1[0],self.coo1[1],self.coo1[0]-self.arrow_F_D_form[0],self.coo1[1]-self.arrow_F_D_form[1])

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
        
        canvas_1.itemconfigure(self.text_vel,text="m = "+str(round(self.target.m,2))+" kg"+"\nI = "+str(round(self.target.I,2))+" kg/m^2"+"\nx, y = "+str(round(self.coo0[0]))+", "+str(round(self.coo0[1]))+"\ntheta_1 = "+str(round(180/pi*self.theta1))+" degrees"+"\nvel = "+str(round(self.absvel,1))+" m/s"+"\nomega = "+str(round(self.omega,3))+" rad/s"+"\nvel_t = "+str(round(self.omega*self.target.r,1))+" m/s"+"\nE_k = "+str(round(self.Ektrans+self.Ekrot))+" J")
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
        canvas_1.itemconfigure(self.text_distint,text=round(self.distint))

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

obj.append(Object_Ball(canvas_1,480,90,15,rho_rubber,0.5,0.6))
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

obj.append(Object_FixedLine(canvas_1,320,200,415,260,5,0.5,0.6))
obj.append(Object_FixedLine(canvas_1,500,550,630,500,5,0.5,0.06))

obj.append(Object_FixedLine(canvas_1,315,h_cnv-100,800,h_cnv-100,5,0.5,0.3))

obj.append(Object_Trace(canvas_1,obj[5],10,20))

#mouse
obj.append(Object_FixedPoint(canvas_1,400,400))
obj.append(Object_Ball(canvas_1,400,450,14,rho_rubber,0.5,0.6))

obj.append(Object_Polygon(canvas_1,550,450,[-35,-15,-12,-27,35,-15,35,15,-35,15],40,rho_rubber,0.5,0.5,0.1,0.1)) #16

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
con.append(Constraint_Distance(obj[14],obj[15],2))

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

obj.append(Object_Line(canvas_1,290,150,360,140,10,40,rho_steel,0.25,0.2)) #270,150 | 390,70

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

obj.append(Object_ShowDistance_Point_LineSegment(obj[23],"coo1",obj[12])) # ändra till ändpunkt 1
obj.append(Object_ShowDistance_Point_LineSegment(obj[23],"coo1",obj[12])) # ändra till ändpunkt 2
obj.append(Object_ShowPhysics(canvas_1,obj[23]))

obj.append(Object_Line(canvas_1,400,75,310,95,3.5,50,rho_rubber,0.5,0.5)) #270,150 | 390,70

obj.append(Object_ShowPhysics(canvas_1,obj[36]))

obj.append(Object_Line(canvas_1,350,500,350,600,3.5,50,rho_rubber,0.5,0.5)) #270,150 | 390,70
obj.append(Object_Line(canvas_1,355,350,350,450,3.5,50,rho_rubber,0.5,0.5)) #270,150 | 390,70

obj.append(Object_Polygon(canvas_1,540,400,[-25,-20,25,-20,25,20,5,9,-25,20],40,rho_rubber,0.5,0.5,0.2,-0.3)) #40
obj.append(Object_Polygon(canvas_1,550,330,[-35,-10,35,-10,35,10,-35,10],40,rho_rubber,0.5,0.5,-0.3,0.1)) #41

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
# polygon och line går ibland genom varandra utan kontakt
# polygon och ine fastnar ibland i varandra
# gör om alla bounding boxes till fyrkanter istället för cirklar (sqrt är långsam)
# initial angle 
# initial angular velocity
# initial velocity
# alla krafter osv funkar för poly2poly, men antipenetreringsisärflyttningen skickar ibland in polygonerna i varandra. något med normalerna?
# ball to box - ändra så att du räknar ut avstånden till alla ytorna och ta det kortaste, istället för y=kx+m?
# ball to box - boll mot kant går igenom varandra aningen. ej fixat pga lathet. se om du orkar. ClosestPointOnLineSegmentPerpDist.
# Det är något skumt med box - line kollisionerna. Ibland går de igenom varandra.
# Zoomfunktion av planen (STORT):
#   - Gör om alla ritkommandon så att de är funktioner med skalning inbyggt som kallas
#   - Lägg till en "måttstock" längst ned på skärmen som skalas med skalningsnivån (tänk ANSYS)
# Ha alla vertices for ett object i x- och y-listor, istället för som separata variabler. Bör ge snyggare kod och förenkla collision- och contact-funktionerna. En enda sådan funktion för alla fall?
# Fixa luftmotståndspilarnas riktning som pekar mest 90 grader isär pga anv av tanh på F_D:s komposanter
# Skapa allmän apply force-funktion som går att använda på vilken punkt som helst på vilket objekt som helst
# Skapa mediemotstånd (luftmotstånd etc.) för rotation
# Varför funkar den nya impulsekvationen i L2FL och inte den gamla?
# Ha kraftpilar osv som ett separat objekt. Som Trace
# cube 
# för fyrkanter/avlånga objekt osv, beräkna dess tvärsnittsarea som möter vinden varje tidssteg - påverkar luftmotsståndet - kan skapa vingar osv
# skapa musdrag av prylar
# skapa rep, med inställning av antalet delar och längd
# allmän kodförbättring - förbättra integrationen så att den är mer samlad?
# skapa portaler?


Timestep()

win_1.mainloop()
