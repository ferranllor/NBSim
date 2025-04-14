a = str(input("Average position: "))[1:-1].split(",")
b = str(input("Body-0  position: "))[1:-1].split(",")

au = [float(a[0]), float(a[1]), float(a[2])]
aok = [0.002313,0.003302,-0.003792]

bu = [float(b[0]), float(b[1]), float(b[2])]
bok = [13.018353,38.377743,23.593529]

print("Var avg: " + str([abs(au[i] - aok[i]) for i in range(3)]))
print("Var b-0: " + str([abs(bu[i] - bok[i]) for i in range(3)]))