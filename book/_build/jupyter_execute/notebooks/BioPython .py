import Bio 

Bio.__version__

from Bio.Seq import Seq

s1 = Seq('ATGGCTTTATTTTCCCGGGA')

s1.complement() 

s1.reverse_complement() 

s1.transcribe() 

s1.back_transcribe() 

s1.back_transcribe() == Seq('ATGGCTTTATTTTCCCGGGA')

s1.translate() 

