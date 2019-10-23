#Sara Jaramillo Cárdenas
#this is the Gene class
#the gene class has to read the gene_information file and keep them, 
#the important thing here is to use the gene_id and the  gene_name as a key for the other Databases
#we have to create a new property call linked but it could be easier to save the info from Cross class


require "./SeedStock.rb"
require "./Cross.rb"

class Gene
    attr_accessor :gene_id
    attr_accessor :gene_name
    attr_accessor :mutant_Genotype
    attr_accessor :linked
    @@Genesa=Hash.new
    
   
    def initialize(params = {})
        @gene_id = params.fetch(:Gene_id,"AT1G69120")
        @gene_name = params.fetch(:Gene_name,"xxx")
        @mutant_Genotype = params.fetch(:Mutant_Genotype,"aa..aa")
        @linked= "" #in this property i want to save the info about gene linked
        @@Genesa[gene_id]=self
       
    end
    def self.guardadata(file_name) #this is for load from file
        
        xf=File.open(file_name,"r")
		primera_linea=xf.readline
		xf.each do |line|
            a,b,c = line.split("\t")
            
            #i wanna check if the gene_id is in the correct format
            regexp = Regexp.new(/A[Tt]\d[Gg]\d\d\d\d\d$/)
            if regexp.match(a) 
                
            else 
                abort "#{a} is not correct"
            end
            Gene.new(:Gene_id =>a, :Gene_name=>b, :Mutant_Genotype =>c)
            
        end

        xf.close

       
    end

      def get_gene
        return @@Genesa
    end
    def Gene.get_gene(var)
        return @@Genesa[var]
    end
   

end

