# This is a Golly Python script for generating the initial rows
# when extending a partial result using zfind, qfind, or gfind.
# Instructions:
# 1. Orient the partial result so that it is traveling upwards.
# 2. Select the row that you wish to extend from (for symmetric
#    searches, only select the left half of the row). The search
#    will essentially use that row and the row below it as the
#    initial rows for the search.
# 3. Run the script and input the period and translation.
# 4. The rows will be saved to "initrows.txt", which can be used
#    directly in a zfind or qfind search.
# 5. To extend with gfind, the rows need to be converted into
#    binary numbers (replace "o" with "1" and "." with "0") and
#    entered into the code according to these instructions:
#    http://www.conwaylife.com/forums/viewtopic.php?f=9&t=925#p6776
#    Note: many compilers allow you to use the "0b" prefix for
#    writing numbers in binary.

from glife import rect, validint
import golly as g

def getrow(x, y, width):
    s = ""
    for i in range(width):
        if g.getcell(x + i, y): s += "o"
        else: s += "."
    return s

firstrow = rect(g.getselrect())
if len(firstrow) == 0: g.exit("There is no selection.")
if firstrow.height != 1: g.exit("Incorrect selection dimensions (height must be 1)")
width = firstrow.width

period = g.getstring("Enter the period.")
if not validint(period): g.exit("Period should be a positive integer.")
period = int(period)
if period <= 0: g.exit("Period should be a positive integer.")

offset = g.getstring("Enter the offset (translation).")
if not validint(offset): g.exit("Offset should be a positive integer.")
offset = int(offset)
if offset <= 0: g.exit("Offset should be a positive integer.")
if period <= offset: g.exit("Period must be greater than offset.")

backOff = [-1 for i in range(period)]

i = 0
while True:
    j = offset
    while ((backOff[(i+j)%period] >= 0) and (j < period)): j += 1
    if j == period:
        backOff[i] = period - i
        break
    backOff[i] = j
    i = (i+j)%period

rows = ["" for i in range(2*period)]
k = 0
mp = 0
x = firstrow.x
y = firstrow.y
for i in range(period):
    rows[k] = getrow(x,y - mp,width)
    rows[k + period] = getrow(x,y + 1 - mp,width)
    k += backOff[k]
    if k >= period:
        k -= period
        mp += 1
    g.run(1)

try:
    f = open('initrows.txt', 'w')
    for i in range(2 * period): f.write(rows[i] + "\n")
    f.close()
except:
    g.warn("Unable to save file.")