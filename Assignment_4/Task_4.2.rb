#Sara Jaramillo
#Bonusss
#i would search the gene funtions and compare between the Best reciprocal hits answers

#'./st.rb'

require  'bio'

require 'stringio'

file_1="TAIR10".to_s
file_2="p.fa".to_s

def create_orthologs_file(filename)
  if File.exists?(filename)
    File.delete(filename)
  end
  return (File.open(filename, "w"))
end
new_file=create_orthologs_file("orthologs_2.txt")
new_file.puts "orthologue pairs between species Arabidopsis and S. pombe"

#create the DataBase Folder
#first delete the older one
#next create a new folder
system("rm -r Databases")
system("mkdir Databases")
system("makeblastdb -in '#{file_1}' -dbtype 'nucl' -out ./Databases/#{file_1}")
system("makeblastdb -in '#{file_2}' -dbtype 'prot' -out ./Databases/#{file_2}")

#Blast
#tblastn --> nucleotides
#blastx --> proteins
E_value = 10**-6
#Ward, N., & Moreno-hagelsieb, G. (2014). Quickly Finding Orthologs as Reciprocal Best Hits with BLAT , LAST , and UBLAST : How Much Do We Miss ?,
# 9(7), 1–6. https://doi.org/10.1371/journal.pone.0101850




factory_f1 = Bio::Blast.local('tblastn', "./Databases/#{file_1}")
factory_f2 = Bio::Blast.local('blastx',"./Databases/#{file_2}")

#fastaformat for reading sequence per sequence
f1_fasta = Bio::FastaFormat.open(file_1)
f2_fasta = Bio::FastaFormat.open(file_2)
fasta_seqs=Hash.new
f1_fasta.each do |seq_f1|

  fasta_seqs[(seq_f1.entry_id).to_s] = (seq_f1.seq).to_s
end
quantity=0

f2_fasta.each do |search_fas_2| # We iterate over each sequence in the search_file
  report = factory_f1.query(search_fas_2)
  if report.hits[0]
    if report.hits[0].evalue <= E_value
      search_id = (search_fas_2.entry_id).delete("\n").to_s

      f1_id= report.hits[0].definition.split("|")[0].delete(" ").to_s #obtein the id

      query_f1=fasta_seqs[f1_id]

      report2 = factory_f2.query(query_f1)
      if report2.hits[0]
        f2_id= report2.hits[0].definition.split("|")[0].delete(" ").to_s

        if search_id == f2_id
       
         
            new_file.puts "#{f1_id}\t-->\t#{f2_id}"
            quantity+=1

         
        end
      end
    end
  end
end

new_file.puts "there are #{quantity} orthologues"



