number = re.split('M|D|I|S', a)
number.pop()
tag = re.findall('[MDIS]', a)
number = [int(x) for x in number]



cumulat = 0
for i in range(len(number)):
  if tag[i] == "D":
    pass
  else:
    cumulat = cumulat + number[i]
    if number[i] > 600:
      length = number[i]
      seq = b[cumulat:(cumulat+length)]
      names = tag[i] + "-" + str(cumulat) + ":" + str(number[i]) + ":"
      print names
      print seq


