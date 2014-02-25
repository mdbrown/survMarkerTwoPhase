
import sys
import commands 


def main():
  
  for batch in range(10):
    batch = batch + 1
    cmd1 = 'sbatch --cpus-per-task=1 --time=3-0 --wrap=\"R --no-save --no-restore --args batch=%d < runNCC_NP.R\"' % (batch)
    (status, output) = commands.getstatusoutput(cmd1)
    print(cmd1)
      
    if status:    ## Error case, print the command's output to stderr and exit
      sys.stderr.write(output)
      sys.exit(1)      
    print output  ## Otherwise do something with the command's output 
 
if __name__ == "__main__":
  main()
