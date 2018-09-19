from Bio import Entrez
from Bio import SeqIO
Entrez.email = 'jinghuawu@berkeley.edu'
orgs = ['Drosophila', 'E.coli', 'human']

glycolysis = ['pyruvate kinase', 'enolase', 
              'Phosphoglycerate mutase', 'phosphofructokinase']

TCA = ['malate dehydrogenase', 'citrate synthase', 
       'aconitase', 'isocitrate dehydrogenase']

pentose_phosphate = ['glucose 6-phosphate dehydrogenase', 'ribose 5-phosphate isomerase', 
                     'transketolase', 'transaldolase']

enzymes = glycolysis + TCA + pentose_phosphate


#gene data
organism = []
ids = []
#get ids
for o in orgs: #for each organism
	for e in enzymes:
		handle = Entrez.esearch(db='nucleotide', #can be anything on Entrez (nucleotide, gene, protein, genome)
								term= o + '[Orgn]' + ' AND ' + e + '[Prot]', #search the term, browser format
								sort='relevance',
								idtype='acc', #types of record IDs are returned (acc=accession number) 
								retmax=1)
		for i in Entrez.read(handle)['IdList']:
			ids.append(i)
			organism.append(o)


description = []
ntseq = []
for i in ids:			
	handle = Entrez.efetch(db='nucleotide', 
                           id=i, #the first one is the most relevant one
                           rettype='gb', #Retrieval type. 
                           retmode='text')
	record = SeqIO.read(handle, "gb")
	description.append(str(record.description))
	ntseq.append(str(record.seq[:30]))

gene_tuple = []
for n in range(len(ids)):
	gene_tuple.append(tuple([ids[n], description[n], organism[n], ntseq[n]]))

#id_tuple = [tuple([i[:]]) for i in ids]
#organism_tuple = [tuple([o[:]]) for o in organism]
#description_tuple = [tuple([d[:]]) for d in description]
#ntseq_tuple = [tuple([nt[:]]) for nt in ntseq]
#gene_tuple = [tuple(g[:]) for g in gene_tuple]

####
import sqlite3 #provide interface
conn = sqlite3.connect('my.db') #create a Connection object that represents the database. ('example.db')
c = conn.cursor() #create cursor object for method calls later.

#gene table
c.execute(
    """CREATE TABLE gene (id INT,
                          name TEXT,
                          description TEXT,
                          organism TEXT,
                          ntseq VARCHAR(20));""")
for i in id_tuple:
	temp = i
	c.execute("""INSERT INTO gene (id) VALUES (?);""", temp)
	conn.commit()

for o in organism_tuple:
	temp = o
	c.execute("""INSERT INTO gene (organism) VALUES (?);""", temp)
	conn.commit()

for d in description_tuple:
	temp = d
	c.execute("""INSERT INTO gene (organism) VALUES (?);""", temp)
	conn.commit()

for nt in ntseq_tuple:
	temp = nt
	c.execute("""INSERT INTO gene (organism) VALUES (?);""", temp)
	conn.commit()


#pathway table
c.execute("""CREATE TABLE pathway (name TEXT, description TEXT)""")
c.execute("""INSERT INTO pathway (name, description) VALUES ('glycolysis', 'a metabolic process that occurs during aerobic and anaerobic respiration of living organisms within the cytoplasm.');""")
c.execute("""INSERT INTO pathway (name, description) VALUES ('TCA', 'a metabolic process that occurs during aerobic and anaerobic respiration of living organisms within the cytoplasm.');""")
c.execute("""INSERT INTO pathway (name, description) VALUES ('pentose_phosphate_pathway', 'hexose monophosphate shunt) is a metabolic pathway parallel to glycolysis.');""")
conn.commit()

#enzyme table
c.execute("""CREATE TABLE enzyme (name TEXT, function TEXT, enzyme_commission TEXT)""")
c.execute("""INSERT INTO enzyme (name, function, enzyme_commission) 
			 VALUES  ('pyruvate kinase', 'catalyzes the final step of glycolysis', '2.7.1.40', 'glycolysis'),
			 		 ('enolase', 'metalloenzyme responsible for the catalysis of the conversion of 2-phosphoglycerate (2-PG) to phosphoenolpyruvate (PEP)', '4.2.1.11', 'glycolysis'),
			 		 ('Phosphoglycerate mutase', 'any enzyme that catalyzes step 8 of glycolysis', '5.4.2.11', 'glycolysis'),
			 		 ('phosphofructokinase', 'a kinase enzyme that phosphorylates fructose 6-phosphate in glycolysis', '2.7.1.11', 'glycolysis'),
			 		 ('malate dehydrogenase', 'an enzyme that reversibly catalyzes the oxidation of malate to oxaloacetate', '1.1.1.37', 'TCA'),
			 		 ('citrate synthase', 'pace-making enzyme in the first step of the citric acid cycle', '2.3.3.1', 'TCA'), 
			 		 ('aconitase', 'an enzyme that catalyses the stereo-specific isomerization of citrate to isocitrate via cis-aconitate in the tricarboxylic acid cycle', '4.2.1.3', 'TCA'), 
			 		 ('isocitrate dehydrogenase', 'an enzyme that catalyzes the oxidative decarboxylation of isocitrate', '1.1.1.42', 'TCA'),
			 		 ('glucose 6-phosphate dehydrogenase', 'a cytosolic enzyme that catalyzes D-glucose 6-phosphate', '1.1.1.49', 'pentose_phosphate_pathway'),
			 		 ('ribose 5-phosphate isomerase', 'catalyzes the conversion between ribose-5-phosphate (R5P) and ribulose-5-phosphate (Ru5P)', '5.3.1.6', 'pentose_phosphate_pathway'),
                     ('transketolase', 'catalyzes two important reactions, which operate in opposite directions in these two pathways', '2.2.1.1', 'pentose_phosphate_pathway'),
                     ('transaldolase', 'an enzyme (EC 2.2.1.2) of the non-oxidative phase of the pentose phosphate pathway', '2.2.1.2', 'pentose_phosphate_pathway');""")
conn.commit()
c.execute("SELECT * FROM enzyme;")
print(c.fetchall())

####
c.execute("SELECT * FROM pathway;")
print(c.fetchone())


