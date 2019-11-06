#class Gene for the Gene id and Int Act id

require 'rest-client'
require 'json'

class Gene
    attr_accessor :gene_id
    attr_accessor :intact_id
    @@GeneHash=Hash.new
    @@GVector=[]

    def initialize(params={})
        @gene_id = params.fetch(:gene_id,"ATddxxxxx")
        @intact_id=params.fetch(:intact_id,nil)
        @@GeneHash[@gene_id]=self
    end

    def self.get_intact_from_togo(x) #this is for get the intact code from togo
        address = "togows.org/entry/ebi-uniprot/#{x}/dr.json"
        response = RestClient::Request.execute(method: :get, url: address)   
        data = JSON.parse(response.body)
        if data[0]['IntAct']
            return data[0]['IntAct'][0][0] 
        else
            return nil
        end
    end

   
        


    end
    
    def self.Open_file(name_file) #open the file
        
        f = File.open(name_file, "r")
        f.each_line do |line|
            line.delete!("\n")
            line =line.sub(/^AT/, "At").to_s #change T -->t
            intactid=get_intact_from_togo(line)
            @@GVector << line #keep all the vlines in a vector
            
            Gene.new(:gene_id=>line,:intact_id=>intactid)
            
        end
        f.close
    end 

    def self.gene_hash
        return @@GeneHash
    end
    def Gene.all_genes_id
                    
     return @@GVector
    end
        
end
Gene.Open_file("ArabidopsisSubNetwork_GeneList.txt")
#puts Gene.all_genes_id

