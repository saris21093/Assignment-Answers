#Sara Jaramillo Cárdenas
#this is the Cross Class
# this class has to do the chi-square
# this class has to identify who is genetically linked 
# it need to take some properties of the other class
# in this class linked has to be filled


require '.\Gene'
require '.\SeedStock'

class Cross
    attr_accessor :parent1
    attr_accessor :parent2
    attr_accessor :f2_wild
    attr_accessor :f2_p1
    attr_accessor :f2_p2
    attr_accessor :f2_p1p2
    @@almacen=Hash.new #i want to keep my info here

    def initialize ( params = {})
        @parent1 = params.fetch(:p1,"xxx")
        @parent2 = params.fetch(:p2,"xxx")
        @f2_wild = params.fetch(:f2wild,"xxx")
        @f2_p1 = params.fetch(:f2p1,"xxx")
        @f2_p2 = params.fetch(:f2p2,"xxx")
        @f2_p1p2 = params.fetch(:f2p1p2,"xxx")
        

        @@almacen[parent1]=self 
    end

    def Cross.readfile(file_name)
        xf=File.open(file_name,"r")
		primera_linea=xf.readline #i took the first line
		#aqui van las otras lineas
		xf.each do |line|
            a,b,c,d,e,f = line.split("\t")
            
            #creating instances 
            Cross.new(:p1 =>a, :p2=>b, :f2wild =>c.to_i, :f2p1 =>d.to_i, :f2p2 =>e.to_i, :f2p1p2 =>f.to_i)
            
        end
        xf.close
    end

    def Cross.chi_square 
        
        @@almacen.each_value do |ins|
    #---------------------------
            a1=ins.f2_wild
            a2=ins.f2_p1
            a3=ins.f2_p2
            a4=ins.f2_p1p2
            total=a1+a2+a3+a4
            y1=(total*0.562)
            y2=(total*0.187)
            y3=(total*0.187)
            y4=(total*0.062)

            r1=((a1-y1)**2)/y1
            r2=((a2-y2)**2)/y2
            r3=((a3-y3)**2)/y3
            r4=((a4-y4)**2)/y4
            sumap=r1+r2+r3+r4
            #puts sumap.round(3)
    #-----------------------------
            p005=7.8147 #this value corresponds to degrees of freedom = 3
            if sumap>=7.8147
                
                #these lines are for save the info in the property of the class Gene call linked
                Gene.get_gene(SeedStock.get_seed_stock(ins.parent1).gene_id).linked = "#{Gene.get_gene(SeedStock.get_seed_stock(ins.parent1).gene_id).gene_name} is genetically linked to #{Gene.get_gene(SeedStock.get_seed_stock(ins.parent2).gene_id).gene_name}"
                Gene.get_gene(SeedStock.get_seed_stock(ins.parent2).gene_id).linked = "#{Gene.get_gene(SeedStock.get_seed_stock(ins.parent2).gene_id).gene_name} is genetically linked to #{Gene.get_gene(SeedStock.get_seed_stock(ins.parent1).gene_id).gene_name}"
                
                
                puts "#{Gene.get_gene(SeedStock.get_seed_stock(ins.parent1).gene_id).gene_name} is genetically linked to #{Gene.get_gene(SeedStock.get_seed_stock(ins.parent2).gene_id).gene_name} with a chi-square score #{sumap}"
            end
        end
        
    end

    def getalmacen
        return @@almacen
    end
    def Cross.getalmacen(var)
        return @@almacen[var]
    end
end
