with open('huh') as f:
	content = f.readlines()
content = [x.strip() for x in content] 
for i in range(2000):
	for j in content:
		if int(i) == int(j):
			print(i+1)
			break
