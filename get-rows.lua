-- This is a Golly Lua script for generating the initial rows
-- when extending a partial result using zfind, qfind, or gfind.
-- Instructions:
-- 1. Orient the partial result so that it is traveling upwards.
-- 2. Select the row that you wish to extend from (for symmetric
--    searches, only select the left half of the row). The search
--    will essentially use that row and the row below it as the
--    initial rows for the search.
-- 3. Run the script and input the period and translation.
-- 4. The rows will be saved to "initrows.txt", which can be used
--    directly in a zfind or qfind search.
-- 5. To extend with gfind, the rows need to be converted into
--    binary numbers (replace "o" with "1" and "." with "0") and
--    entered into the code according to these instructions:
--    https://www.conwaylife.com/forums/viewtopic.php?p=6776#p6776
--    Note: many compilers allow you to use the "0b" prefix for
--    writing numbers in binary.

local g = golly()
local i

local function getrow(x, y, width)
   local s = ""
   for i = 0, width-1, 1 do
      if g.getcell(x + i, y) == 1 then
         s = s .. "o"
      else
         s = s .. "."
      end
   end
   return s
end

local firstRow = g.getselrect()
if firstRow == nil then
   g.exit("There is no selection")
elseif firstRow[4] ~= 1 then
   g.exit("Incorrect selection dimensions (height must be 1)")
end
local width = firstRow[3]

local period = tonumber(g.getstring("Enter the period."))
if period == nil or period ~= math.floor(period) or period <= 0 then
   g.exit("Period should be a positive integer.")
end

local offset = tonumber(g.getstring("Enter the offset (translation)."))
if offset == nil or offset ~= math.floor(offset) or offset <= 0 then
   g.exit("Offset should be a positive integer.")
end

if offset >= period then
   g.exit("Offset should be less than period.")
end

local backOff = {}
for i = 0, period-1, 1 do
   backOff[i] = -1;
end

i = 0
while true do
   local j = offset
   while ((backOff[(i+j)%period] >= 0) and (j < period)) do
      j = j + 1
   end
   if j == period then
      backOff[i] = period - i
      break
   end
   backOff[i] = j
   i = (i+j)%period
end

local k = 0
local mp = 0
local x = firstRow[1]
local y = firstRow[2]
local rows = {}

for i = 1, period, 1 do
   rows[k] = getrow(x,y - mp,width)
   rows[k + period] = getrow(x,y + 1 - mp,width)
   k = k + backOff[k]
   if k >= period then
      k = k - period
      mp = mp + 1
   end
   g.run(1)
end

local f = io.open("initrows.txt", "w")
if f then
   for i = 0, 2*period-1, 1 do
      f:write(rows[i] .. "\n")
   end
   f:close()
else
   g.warn("Unable to save file.")
end