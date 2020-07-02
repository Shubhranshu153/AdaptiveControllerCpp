import matplotlib.pyplot as plt
import csv
import numpy as np
import matplotlib.animation as animation
from matplotlib.patches import Circle 
import pdb




if __name__ == '__main__':
    
    data={'time':[], 'mass1':[],'mass2':[],'q1':[],'q2':[],'q1dot':[],'q2dot':[],'Dq1':[],'Dq2':[],'Dq1dot':[],'Dq2dot':[]}
    with open('q1.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            data['time'].append(float(row[0]))
            data['q1'].append(float(row[1]))
            data['q2'].append(float(row[2]))
            data['q1dot'].append(float(row[3]))
            data['q2dot'].append(float(row[4]))

            data['mass1'].append(float(row[5]))
            data['mass2'].append(float(row[6]))

            data['Dq1'].append(float(row[7]))
            data['Dq2'].append(float(row[8]))
            data['Dq1dot'].append(float(row[9]))
            data['Dq2dot'].append(float(row[10]))
    
    fig1=plt.figure(1)
    plt.title('q1 vs q1 desired')
    plt.plot(data['time'],data['q1'], label='q1',alpha=0.8)
    plt.plot(data['time'],data['Dq1'],'--r',label='desired q1')
    plt.legend(loc="upper left")

    fig2=plt.figure(2)
    plt.title('q2 vs q2 desired')
    plt.plot(data['time'],data['q2'],label='q2',alpha=0.8)
    plt.plot(data['time'],data['Dq2'],'--r',label='desired q2')
    plt.legend(loc="upper left")
    
    # plt.show()


    def make_plot(i):
        # Plot and save an image of the double pendulum configuration for time
    # point i.
    # The pendulum rods.
        ax.plot([0, x1[i], x2[i]], [0, y1[i], y2[i]], lw=2, c='k')
        # Circles representing the anchor point of rod 1, and bobs 1 and 2.
        c0 = Circle((0, 0), r/2, fc='k', zorder=10)
        c1 = Circle((x1[i], y1[i]), r, fc='b', ec='b', zorder=10)
        c2 = Circle((x2[i], y2[i]), r, fc='r', ec='r', zorder=10)
        ax.add_patch(c0)
        ax.add_patch(c1)
        ax.add_patch(c2)

        # The trail will be divided into ns segments and plotted as a fading line.
        ns = 20
        s = max_trail // ns

        for j in range(ns):
            imin = i - (ns-j)*s
            if imin < 0:
                continue
            imax = imin + s + 1
            # The fading looks better if we square the fractional length along the
            # trail.
            alpha = (j/ns)**2
            ax.plot(x2[imin:imax], y2[imin:imax], c='r', solid_capstyle='butt',
                    lw=2, alpha=alpha)

        # Centre the image on the fixed anchor point, and ensure the axes are equal
        ax.set_xlim(-L1-L2-r, L1+L2+r)
        ax.set_ylim(-L1-L2-r, L1+L2+r)
        ax.set_aspect('equal', adjustable='box')
        plt.axis('off')
        plt.savefig('frames/_img{:04d}.png'.format(i), dpi=72)
        plt.cla()
    

    # pdb.set_trace()
    x1 = 0.5 * np.sin(data['q1'])
    y1 = -0.5 * np.cos(data['q1'])
    x2 = x1 + 1 * np.sin(np.array(data['q1'])+np.array(data['q2']))
    y2 = y1 - 1 * np.cos(np.array(data['q1'])+np.array(data['q2']))
    r = 0.05
    trail_secs = 1
    dt=0.01
    # This corresponds to max_trail time points.
    max_trail = int(trail_secs / dt)
    L1=0.5
    L2=1
 
   

    fps = 10
    di = int(1/fps/dt)
    fig = plt.figure(figsize=(8.3333, 6.25), dpi=72)
    ax = fig.add_subplot(111)
    i=0
    while i<np.array(data['time']).shape[0]:
        
        print(i , '/', np.array(data['time']).shape[0])
        make_plot(i)
        i=i+100
