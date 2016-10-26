import sys

opts = {
  'window': 100000,
  'intrachrom': False
}
args = sys.argv[1:]
for a, arg in enumerate(args):
  if arg.startswith('-'):
    if arg.startswith('-w') or arg.startswith('--window'):
      opts['window'] = int(args[a+1])
    if arg.startswith('-i'):
      opts['intrachrom'] = True
  else:
    args = args[a+1:]
    break

def main():
  breakpoints =[]
  ensemble = []
  for filename in args:
    breakpoints += get_breakpoints(filename)
  ensemble = find_fusions(breakpoints)
  ensemble.sort(key=len, reverse=True)
  print '%s ensemble fusions found:' % len(ensemble)
  for calls in ensemble:
    print ''
    for call in calls:
      result = '%s:%s %s >< %s %s:%s (%s, %s:%s)' % (call['5p']['chr'], call['5p']['pos'], call['5p']['gene'], call['3p']['gene'], call['3p']['chr'], call['3p']['pos'], call['caller'], call['file'], call['l'])
      print result


def find_fusions(breakpoints):
  fusions = []
  for p1, point1 in enumerate(breakpoints):
    match = [point1]
    for point2 in breakpoints[p1+1:]:
      if overlap(point1, point2):
        match.append(point2)
    callers = map(lambda x: x['caller'], match)
    if len(set(callers)) > 1:
      for calls in fusions:
        if set(map(id, match)) < set(map(id, calls)):
          break
      else:
        fusions.append(match)
  return fusions


def get_breakpoints(filename):
  # Returns list of positions, e.g. [ [[1, 123], [2, 789]], [[3, 321], [4, 987]] ]
  caller = 0
  breakpoints = []
  #print 'Retrieving breakpoints for %s' % filename
  for l,line in enumerate(open(filename)):
    if caller:
      point = {'caller': caller, 'line': line, 'file': filename, 'l': l}
      cols = line.split('\t')
      if caller == 'STAR':
        point['5p'] = decolonize(cols[4])
        point['3p'] = decolonize(cols[7])
        point['5p']['gene'], point['3p']['gene'] = cols[0].split('--')
      elif caller == 'FUSION_CATCHER':
        point['5p'] = decolonize(cols[8])
        point['3p'] = decolonize(cols[9])
        point['5p']['gene'] = cols[0]
        point['3p']['gene'] = cols[1]
      elif caller == 'DEFUSE':
        point['5p'] = {'chr': cols[24], 'pos': cols[37], 'gene': cols[30]}
        point['3p'] = {'chr': cols[25], 'pos': cols[38], 'gene': cols[31]}
      elif caller == 'CHIMERA':
        pos1 = cols[2] if cols[8] == '-' else cols[1] # if 5' -> 3', use strand1 end
        pos2 = cols[4] if cols[9] == '-' else cols[5] # if 5' -> 3', use strand2 start
        point['5p'] = {'chr': cols[0].split('chr')[1], 'pos': pos1, 'gene': cols[12]}
        point['3p'] = {'chr': cols[3].split('chr')[1], 'pos': pos2, 'gene': cols[13]}
      point['5p']['gene'] = point['5p']['gene'].upper()
      point['3p']['gene'] = point['3p']['gene'].upper()
      breakpoints.append(point)
    elif line.startswith('#chrom5p'): caller = 'CHIMERA'
    elif line.startswith('#fusion_name'): caller = 'STAR'
    elif line.startswith('Gene_1_symbol'): caller = 'FUSION_CATCHER'
    elif line.startswith('cluster_id'): caller = 'DEFUSE'
  return breakpoints


def id(point):
  return (point['file'], point['l'])


def decolonize(vals):
  # converts "chr4:1234:+" to {'chr': '4', 'pos': '1234'}
  chrom, pos = vals.split(':')[:2]
  if chrom.startswith('chr'):
    chrom = chrom.split('chr')[1]
  return {'chr': chrom, 'pos': pos}


def overlap(pt1, pt2):
  p1c1, p1c2, p2c1, p2c2 = pt1['5p']['chr'], pt1['3p']['chr'], pt2['5p']['chr'], pt2['3p']['chr']
  if p1c1 == p2c1 and p1c2 == p2c2:
    if p1c1 != p1c2 or opts['intrachrom']:
      if abs(int(pt1['5p']['pos']) - int(pt2['5p']['pos'])) <= opts['window'] and abs(int(pt1['3p']['pos']) - int(pt2['3p']['pos'])) <= opts['window']:
        return True
  return False


if __name__ == '__main__': main()
