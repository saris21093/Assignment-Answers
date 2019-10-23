#Sara Jaramillo Cárdenas
#Task1
#this is the SeedStock Class


class SeedStock
	
	attr_accessor :gene_id
	attr_accessor :seed_stock_id
	attr_accessor :last_planted
	attr_accessor :grams_remaining
    attr_accessor :storage
    @@Dicc=Hash.new #i want to keep all my properties here


	def initialize (params = {}) #voy a meter aqui un hash con todos los parametros
		#como seed_stock, gene_id Last_planted Storage y Grams reimaining
		@seed_stock_id = params.fetch(:Seed_Stock, "Xxxx")
		@gene_id = params.fetch(:Mutant_Gene_ID,"XXxXxxxxx")
		@last_planted = params.fetch(:Last_Planted,"01/01/2000") 
		@storage=params.fetch(:Storage,"XXXXx")
        @grams_remaining = params.fetch(:Grams_Remaining, 0)#it has to be an int
        @@Dicc[seed_stock_id]=self

		
    end
    def SeedStock.load_from_file(file_name)

		xf=File.open(file_name,"r")
		$primera_linea=xf.readline
		#for reading the next lines
		xf.each do |line|
			ss,mgid,lp,s,gr = line.split("\t")
            #i wanna check if the gene_id is in the correct format
            regexp = Regexp.new(/A[Tt]\d[Gg]\d\d\d\d\d$/)
            if regexp.match(mgid) 
            else 
                abort "#{mgid} is not correct"
            end

           #here im creating the instances
            SeedStock.new(:Seed_Stock =>ss, :Mutant_Gene_ID=>mgid, :Last_Planted =>lp, :Storage =>s, :Grams_Remaining=>gr.to_i)
        end
        xf.close
    end

    def sembrar(cantidad) #this class is for plant the seeds, cantidad is the number of grams that i need to plant
        
        if @grams_remaining<=cantidad
            puts "WARNING: we have run out of Seed Stock #{@seed_stock_id}"
            #puts "you dont have that quantity of #{@seed_stock_id} you only have #{@grams_remaining}, so we only plant those grams"
            @grams_remaining=0
        else
            @grams_remaining=@grams_remaining-cantidad
        end
    end

    def self.write_database(new_one)
        if File.exists?(new_one) #i dont want to overwrite, i want to create always a new file
            File.delete(new_one)
        end
        t = Time.now
        fecha= t.strftime("%d/%m/%Y") #this is for get the current date

        xo=File.open(new_one,"a+")
        xo.puts $primera_linea #to write the header in the new file
        @@Dicc.each_value do |ins| #this is for writing all the instances with their properties in the new file
            
            xo.puts "#{ins.seed_stock_id}\t#{ins.gene_id}\t#{fecha} \t #{ins.storage}\t #{ins.sembrar(7)}"
        end
    end
    #these are for access to the database from where i want
    def get_seed_stock
        return @@Dicc
    end
    def SeedStock.get_seed_stock(var)
        return @@Dicc[var]
    end

    


end



