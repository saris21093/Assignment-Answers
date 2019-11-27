# @author Sara Jaramillo Cárdenas
# the inputs of the code is the file with the target genes
# the code is going to match cttctt in the exon's sites
# the code is going to generate a its own feauture
# the outputs of the code are:
# 1. a file with the genes do not  have exons with the cttctt
# 2. two gff files, one in gene coordinates and the other with chromosome coordinates


require 'net/http'
require 'bio'

# function for deleting existing files
# @param file_name [string] file's name with the format

def file_exists(file_name)
    if File.exists?(file_name) 
       File.delete(file_name)
    end
end

# create files and writting the headers

file_genes="genes.gff"
file_chromosoma="chr.gff"
file_no_matches = "No_Matches.txt"
file_exists(file_genes)
file_exists(file_chromosoma)
file_exists(file_no_matches)

fo=File.open(file_genes, "a+")
fo.puts "##gff-version 3"
f1=File.open(file_no_matches,"a+")
f1.puts "The genes where there are not a match"
f2=File.open(file_chromosoma,"a+")
f2.puts "##gff-version 3"


# function for getting data from the file
# @param file [String] file's name with the format
# @return gene [array] list of all the genes in the file

def get_data_from_file(file)
    f = File.open(file, "r")
    gene = Array.new
    f.each_line do |line|
      gene << line.delete("\n") # We remove end of line \n 
    end
    return gene
end

# function for get the data from ebi
# @param gene [array] list of genes
# @return bioseq [Bio::Sequence] the sequence of the gene
# @return chr_id [string] the chromosome's id
# @return chr_cord [array] the chromosome's coordinates

def obtain_data_from_ebi(gene)
    address = URI("http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensemblgenomesgene&format=embl&id=#{gene}") 
    response = Net::HTTP.get_response(address)  # use the Net::HTTP object "get_response" method
    record = response.body
    entry = Bio::EMBL.new(record) 
    ac_seq=  entry.accession.split(":")
    chr_id=ac_seq[2]
    chr_cord=[ac_seq[3].to_i,ac_seq[4].to_i]
    #obtein the secuence
    bioseq = entry.to_biosequence
    return bioseq, chr_id, chr_cord
end
# function for get the targets cttctt in the exons
# @param bioseq_seq [Bio::Sequence] the sequence of the gene
# @return targets [hash] the hash'keys are the coordinates of the matches and the values are the exon_id and the string
def obtain_target_from_seq(bioseq_seq)
    #we are going to keep all the targets from each gene,
    #hash of hashes :)
    targets=Hash.new
    len_bioseq = bioseq_seq.length + 1
    nstrand_targets = []
    pstrand_targets = []
    for i in (0..bioseq_seq.length-5)
        if bioseq_seq.complement[i..i+5] == "cttctt"
            nstrand_targets << [i,i+5]
        end
        if bioseq_seq[i..i+5] == "cttctt"
            pstrand_targets << [i,i+5]
        end
    end
    #im going to obtein the exon sites
    bioseq_seq.features.each do |feature|
        position = feature.position
        #not remote entries
        next unless (feature.feature == 'exon' ) && (not position =~ /[A-Z]/)
        exon_id = feature.qualifiers[0].value.gsub('exon_id=', '')
        #we can find exons in the the complement or the postive strand, so the way for obtaining the data is different
        #i wanna keep the strand forward or reverse
        if position =~/complement/
            exon_site_negative=position.tr('complement()',"")
            aux=exon_site_negative =~ /\./ 
            #we have to change the start and final point, the ebi documentation shows  --> x..y, y is the start point and x is the final point
            #location 1 --> l1 
            #location 2 --> l2
            l1= len_bioseq - exon_site_negative[0,aux].to_i
            l2= len_bioseq - exon_site_negative[aux+2,exon_site_negative.length].to_i
            exon_site_negative=[l2,l1] # they are the location of the exons
            #we are going to check if the target that we did before is in the exon site
            is_inside_exon = check_target_in_exon(exon_id,nstrand_targets,'reverse',len_bioseq,exon_site_negative)
            #is_inside_exon will be a hash with positions, exon_id and the strand for the gff file
            unless is_inside_exon.nil? #if the hash is not nil i could keep it in other hash
                targets = targets.merge(is_inside_exon)
            end
        else
            exon_site_positive=position
            aux=exon_site_positive =~ /\./ 
            l1= exon_site_positive[0,aux].to_i
            l2= exon_site_positive[aux+2,exon_site_positive.length].to_i
            exon_site_positive=[l1,l2]
            is_inside_exon = check_target_in_exon(exon_id,pstrand_targets,'forward',len_bioseq,exon_site_positive)
            #is_inside_exon will be a hash with positions, exon_id and the strand for the gff file
            unless is_inside_exon.nil? #if the hash is not nil i could keep it in other hash
                targets = targets.merge(is_inside_exon)
            end
        end
    end
    return targets
end

# function for check if the target is inside the exon's sites
# @param exon_id [string] the exon's id
# @param strand [string] the direction of the strand
# @param len_bioseq [integer] the length of the sequence
# @param exones_site [array] the coordinates of the exon's site
# @return target_in_exon [hash] the hash'keys are the coordinates of the matches and the values are the exon_id and the string


def check_target_in_exon(exon_id,strand_target,strand,len_bioseq,exones_site)
    #we are going to keep the exon_id, the strand and the location of the target for the gff file
    target_in_exon = Hash.new
    strand_target.each do |pos|
        #c=[]--> vector auxiliar 
        c=[]
        pos.zip(exones_site).map { |x, y| c << x-y}
        if c[0] >=0 && c[1] <= 0
            if strand == 'reverse'
                #the target is in the exon
                #for the format of ebi we have to change again our start and final points 
                e_start = len_bioseq - pos[1].to_i
                e_final = len_bioseq - pos[0].to_i
                target_in_exon[[e_start,e_final]] = [exon_id,'-']
            else
                target_in_exon[[pos[0],pos[1]]] = [exon_id,'+']
            end
        end
    end
    if not target_in_exon.empty? # We check there are targets inside the exon
        return target_in_exon
    end
end

# function for creating new features to the EnsEMBL Sequence object
# @param bioseq [Bio::Sequence] the sequence of the gene
# @param targets [hash] the hash'keys are the coordinates of the matches and the values are the exon_id and the string


def create_features_ensembl_seq_obj(bioseq,targets)
    targets.each do |key,value|
        f1 = Bio::Feature.new("target_CTTCTT","#{key[0]}..#{key[1]}")
        #im no sure if it is a interior coding exon 
        f1.append(Bio::Feature::Qualifier.new('interior coding exon', "#{value[0]}"))
        f1.append(Bio::Feature::Qualifier.new('strand', "#{value[1]}"))
        bioseq.features << f1
    end
end
# function for changing the gene coordinates to chromosome coordinates
# @param target_hash [string] the hash'keys are the coordinates of the matches and the values are the exon_id and the strand
# @param chr_cord [array] the chromosome's coordinates
# @return chr_target_hash [hash] the hash'keys are the coordinates of the matches and the values are the exon_id and the strand

def change_cord(target_hash,chr_cord)
    chr_target_hash ={}
    target_hash.each do |key, value|
        v_cor=[chr_cord[0] + key[0],chr_cord[0] + key[1]]
        chr_target_hash[v_cor]=value
    end
    return chr_target_hash
end

# function for writing in the txt and gff files
# @param gene [array] gene's list
# @param fo [string] gene's gff file
# @param f1 [string] NoMatches's txt file
# @param f2 [string] chromosome's gff file

def write_gff_files(gene,fo,f1,f2)
    bioseq_seq, chr_id, chr_cord = obtain_data_from_ebi(gene) 
    target_hash = obtain_target_from_seq(bioseq_seq)
    create_features_ensembl_seq_obj(bioseq_seq,target_hash)
    chr_target_hash=change_cord(target_hash,chr_cord)
    if target_hash.empty?
        f1.puts "#{gene} \n"
    else
        #for chr gff file, #this is for the parent gene chr_id,. --> source, "gene", the coordinates in the chr, . --> score, strand, . --> phase, gene_id
        f2.puts "#{chr_id}\t.\tgene\t#{chr_cord[0]}\t#{chr_cord[1]}\t.\t+\t.\tID=#{gene}"
    end
    #this is for the chromosomas, chr_id, . --> source, featuretype, cordinates of the target, . --> score, strand, . --> phase, exon id and the Parent identifiers, ID=exon00001;Parent=mrna0001
    chr_target_hash.each do |key,value|
        f2.puts "#{chr_id}\t.\tinterior coding exon\t#{key[0]}\t#{key[1]}\t.\t#{value[1]}\t.\t#{value[0]};Parent=#{gene}"
    end
    # each loop for write in the gff file
    bioseq_seq.features.each do |feature|
        featuretype = feature.feature
        next unless featuretype == "target_CTTCTT"
        position = feature.position
        qual = feature.assoc          
        positionss= position.split("..")
        fo.puts"#{gene}\t.\t#{featuretype}\t#{positionss[0]}\t#{positionss[1]}\t.\t#{qual["strand"]}\t.\t#{qual["interior coding exon"]}"
    end
end
#----------------------------------------

# get all the genes from the file
gene_id = get_data_from_file("ArabidopsisSubNetwork_GeneList.txt")
# each loop for write in the txt and gff file 
gene_id.each do |gene|
    write_gff_files(gene,fo,f1,f2)
  
end