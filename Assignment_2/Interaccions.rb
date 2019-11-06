#Sara Jaramillo 
#this is the class interations 
#This class do everything :(

require 'rest-client'
require 'json'  # to handle JSON format
require './Gene.rb'

class Interaction
    attr_accessor :network
    attr_accessor :members
    attr_accessor :nodos
    attr_accessor :gen_int_lis
    attr_accessor :go
    attr_accessor :kegg
    attr_accessor :keggpath
    attr_accessor :goterm

    @@InteractionHash=Hash.new
    @@Deep=1
    @@variable=[]
    @@Vec=Array.new(@@Deep+1) { Array.new() }
    @@auxd=0
    

    

    def initialize(params ={})
        @network=params.fetch(:network,0)
        @members=params.fetch(:members,[])
        @nodos=params.fetch(:nodos,0)
        @gen_int_lis=params.fetch(:gen_int_lis,[])
        @go=params.fetch(:go, [])
        @goterm=params.fetch(:goterm, [])
        @kegg=params.fetch(:kegg, [])
        @keggpath=params.fetch(:keggpath, [])
        @@InteractionHash[@network]=self
    end

    #---FETCH----
    def self.fetch(url, headers = {accept: "*/*"}, user = "", pass="")
        response = RestClient::Request.execute({
              method: :get,
              url: url.to_s,
              user: user,
              password: pass,
              headers: headers})
            return response
            
            rescue RestClient::ExceptionWithResponse => e
              $stderr.puts e.response
              response = false
              return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
            rescue RestClient::Exception => e
              $stderr.puts e.response
              response = false
              return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
            rescue Exception => e
              $stderr.puts e
              response = false
              return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
    end 
    #---END FETCH---
   
    
    def self.interactors(gene_id) #class method for searching the interaction
        aux=[]
        res = fetch ("http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/interactor/#{gene_id}?species:3702?format=tab25")
        #--------------res--------------------------------------------
        if res
            body = res.body
            lines = body.split("\n")
            #-------iteration of lines from ebi----------------
            lines.each do |i|
                #-------this is for finding locus of p1
                
                a= i.split("\t")[4]
                num_expr = (/A[Tt]\d[Gg]\d\d\d\d\d/) 
                a1 = num_expr.match(a).to_s
                if a1==""
                    next
                else
                    g1 = a1
                end
                #---------end locus =g1----------------
                 #-------this is for finding locus of p2--------
                b=i.split("\t")[5]
                num_expr = (/A[Tt]\d[Gg]\d\d\d\d\d/) # expresión regular 
                b1 = num_expr.match(b).to_s  
                if b1==""
                    next
                else
                    g2 = b1
                end
                #---------end locus =g2----------------------------

                #--------Finding miScore--------------------------
                c= i.split("\t")[14]
                intact_miscore =c.sub(/intact-miscore:/, "").to_f
                
                #--------end miScore = intact_miscore-------------

                #-------find the gene----------------------------
                if gene_id == g1
                    gene=g2
                else
                    gene=g1
                end
                #--------find the gene =gene---------------------
                
                #----------Filters--------------------------------
                case

                when intact_miscore < 0.485 #filter 2: when the miscore is lower than 0.485, this value was taken from the literature
                    next    
                when g1 == g2 #filter 1: when the proteins are the same
                    next 
                when aux.include?(gene) #filter 4: when the protein is already save
                   next
                end
                #-------end of filters----------------------------------
                
                
                cont=0
                unless @@Vec[@@auxd] == nil
                     #puts "aqui es por que es nill entro "
                    @@Vec.each do |k|
                        if k.include?(gene) == FALSE
                            cont+=1
                         #puts "aqui entro por que el gen ya esta en la variable"
                         #puts k
                        end
                         
                    end
                        if cont == @@Vec.length
                            aux << gene
                           
                        end
                    end
                 end#----------end of iteration of lines from ebi
         return aux
        end#------------------end from res----------------------------
    end #end of interactors

    def self.Get_interactors # with a while bucle i will interatc over all the branches that i need
        network=0
        Gene.gene_hash.each_value do |clave|
            unless clave.intact_id == nil
                @@auxd=0
                @@Vec=Array.new(@@Deep+1) { Array.new() }
                @@variable=[]
                h=[]
                h << clave.gene_id
                @@variable << h
                @@Vec[@@auxd] = h
                nodos=1
                a=Interaction.interactors(clave.gene_id)
                unless a==[]
                    @@variable << a
                    @@auxd=1
                    @@Vec[@@auxd] = a
                    c=(a.length)-1
                    nodos+=@@Vec[@@auxd].length
                    while @@auxd <= @@Deep
                        b=Interaction.interactors(@@Vec[@@auxd][c])
                        @@Vec[@@auxd+1]=b
                        @@variable<<b
                        unless c==0
                            c-=1
                        else
                            @@auxd+=1
                            if (@@Vec[@@auxd].length) == 0
                                break
                            else
                            c=(@@Vec[@@auxd].length)-1
                            nodos+=@@Vec[@@auxd].length
                            end
                        end
                    end
                    network+=1
                    var=Interaction.Busqueda
                    go, kegg, keggpath, goterm = Interaction.Get_go_kegg
                    Interaction.new(:network=> network, :nodos=>nodos, :members=>@@Vec, :kegg => kegg, :keggpath => keggpath, :go => go, :goterm => goterm, :gen_int_lis => var.length)


                end
            end
        end

    end

    def self.Busqueda #this search for which genes are include in the list
        var=[]
        @@Vec.each {|i| i.each {|j|   var << j if Gene.all_genes_id.include?(j)}}
        

        return var
    end

    def self.interaction_hash
        return @@InteractionHash
    end

    def self.Get_go_kegg #get the go and the kegg
       
        kegg=[]
        go=[]
        keggpath=[]
        goterm=[]
        @@Vec.each do |k|
            k.each do |j|
                
                addressGO = "http://togows.org/entry/ebi-uniprot/#{j}/dr.json"
                responseGO = RestClient::Request.execute(method: :get, url: addressGO)   
                dataGO = JSON.parse(responseGO.body)
                if dataGO[0]["GO"]
                    dataGO[0]["GO"].each do |i|
                        if i[1] =~ /^P:/ # We must check the GO refers to a biological proccess (it will start with a 'P')
                          goterm << i[1].sub(/P:/, "")
                          go << i[0]
                         end
                    end
                end

                addressKEGG = "http://togows.org/entry/kegg-genes/ath:#{j}/pathways.json"
                responseKEGG = RestClient::Request.execute(method: :get, url: addressKEGG)   
                dataKEGG = JSON.parse(responseKEGG.body)
                if dataKEGG[0]
                    dataKEGG[0].each do |pathid, name_path|
                        kegg << pathid
                        keggpath << name_path
                    end
                end
            end
        end
        return go, kegg, keggpath, goterm
    
            
    end

    def self.write_database(new_one) #writing in a new file
        if File.exists?(new_one) #i dont want to overwrite, i want to create always a new file
            File.delete(new_one)
        end
         xo=File.open(new_one,"a+")
        xo.puts "interacciones" #to write the header in the new file
        @@InteractionHash.each_value do |ins| #this is for writing all the instances with their properties in the new file
            
            xo.puts " Network Number: #{ins.network}. \n Number of Nodos:  #{ins.nodos} \n there are #{ins.gen_int_lis} genes from the list"
            ins.go.each_with_index {|x,y| xo.puts "#{x} \t #{ins.goterm[y]} \n"  }
            ins.kegg.each_with_index {|l,z| xo.puts "#{l} \t #{ins.keggpath[z]} \n"}
            
            
        end
    end

end

Interaction.Get_interactors
Interaction.write_database("new_one.txt")