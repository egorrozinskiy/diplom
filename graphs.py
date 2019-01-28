file= open('tom.txt', 'r')
N= 2049
f = open('toma.txt', 'r')
x = []
y = [0]*N
content = []
number_x = 0
number_y = 1
for line in file:
    x += str(line).split()
print(x)

for i in range(N):
    y[i] = float(x[i])

print(y)




