#Sara Jaramillo Cárdenas
#this is a file which allows us to manage all the DataBases


require '.\Gene'
require '.\SeedStock'
require '.\Cross'

SeedStock.load_from_file("seed_stock_data.tsv")
SeedStock.write_database("new_file.tsv")
Gene.guardadata("gene_information.tsv")
Cross.readfile("cross_data.tsv")
Cross.chi_square
puts "Final Report:"
puts Gene.get_gene("AT5G20240").linked
puts Gene.get_gene("AT1G30950").linked
 
