## Automated bash script that breaks up a prmtop file into smaller 
## PRMTOPS, given the respective residue ranges.
## It uses ANTE-MMPBSA.py, which comes packaged with Antechamber.

## It DOES NOT handle .rst7 file generation, since that is easier to do
## from a larger rst7 file using VMD (for example).

## Usage:
## chainBreaker.sh source_Prmtop_file.prmtop residue_range_01 residue_range_02 ...
## Example Usage:
## bash chainBreaker.sh example.prmtop 1-53 54-104

if [[ $# < 2 ]]; then
    echo "Error: Need to pass at least one residue range for this script to work."
    exit 1
fi

if ! command -v ante-MMPBSA.py &> /dev/null ; then
    echo "Error: ante-MMPBSA.py could not be found. Please install AmberTools (for example using miniconda!)."
    exit 1
fi

fn=`basename $1 .prmtop` 
fileNumber=0
## Iterate over arguments and execute commands for each one
for arg in "${@:2}"; do
    ante-MMPBSA.py -p $1 -r $fn.$fileNumber.prmtop -l temp.prmtop -m \(:$arg\)
    rm temp.prmtop
    fileNumber=$((fileNumber+1))
done
