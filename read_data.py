import matplotlib.pyplot as plt
import csv

if __name__ == '__main__':
    time=[]
    angle=[]
    with open('q1.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            time.append(float(row[0]))
            angle.append(float(row[1]))
        
    plt.plot(time,angle);
    plt.show()
