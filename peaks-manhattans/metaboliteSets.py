import pdb

modelType = 'diet'

if modelType == 'main':
    groupedPeaks_main = open('Z:/Drosophila/Bayesian/general/results/metabolites/updated/plots/groupedPeaks_Main2.csv', 'rb')
    theData = groupedPeaks_main.readlines()

elif modelType == 'diet':
    groupedPeaks_diet = open('Z:/Drosophila/Bayesian/general/results/metabolites/updated/plots/groupedPeaks_Diet2.csv', 'rb')
    theData = groupedPeaks_diet.readlines()

elif modelType == 'comparison':
    groupedPeaks_main = open('Z:/Drosophila/Bayesian/general/results/metabolites/updated/plots/groupedPeaks_Main2.csv', 'rb')
    groupedPeaks_main_data = groupedPeaks_main.readlines()
    groupedPeaks_diet = open('Z:/Drosophila/Bayesian/general/results/metabolites/updated/plots/groupedPeaks_Diet2.csv', 'rb')
    groupedPeaks_diet_data = groupedPeaks_diet.readlines()

    theData = groupedPeaks_main_data + groupedPeaks_diet_data

#metabNetwork_sif = open('Z:/Drosophila/Bayesian/general/results/metabolites/updated/plots/metabNetwork_%s.sif'%modelType, 'w')
#metabNetwork_sif2 = open('Z:/Drosophila/Bayesian/general/results/metabolites/updated/plots/metabNetwork_%s2.tab'%modelType, 'w')
#metabNetwork_sif2.write('Metab1\tMetab1.Type\tInteraction.Type\tMetab2\tMetab2.Type\n')
#metabGMT_file = open('Z:/Drosophila/Bayesian/general/results/metabolites/updated/plots/metabGMT_%s2.gmt'%modelType, 'w')

mainModels = ['Additive', 'Dominant', 'Full', 'Main']
dietModels = ['Add-Diet', 'Dom-Diet', 'Full-Diet', 'Diet']
all_degreeDist = {}

#metabFreq[metabolite] = [itself, withID, withoutID]
metabFreq = {}

for i in xrange(len(theData)):
    if i == 0 or (modelType=='comparison' and i==560):
        #pdb.set_trace()
        continue
    dataSplit = theData[i].split(',')
    
    metabSetDescrip = 'Includes %s metab and %i QTLs with Strength %s' % (dataSplit[10], int(dataSplit[11]), dataSplit[9])
    metabSet = dataSplit[3].split(';')
    model = dataSplit[1][1:-1]
    metabSetName = '%s-%s:%s-%s' % (model, dataSplit[5][1:-1], dataSplit[6], dataSplit[7])
    
    if len(metabSet) == 1:
        theMetab = metabSet[0].replace('"', '').replace('X_', '')
        #metabNetwork_sif.write('%s\n' % theMetab)
        #metabGMT_file.write('%s\t%s\t%s\n' % (metabSetName, metabSetDescrip, metabSet[0][1:-1]))
        try:
            print(int(theMetab))
            metab_type = 'Unidentified'
            #metabNetwork_sif2.write('%s\tUnidentified\t \t \t \n' % theMetab)
        except ValueError:
            print(theMetab)
            if (theMetab[0] == '_'):
                pdb.set_trace()
            metab_type = 'Identified'
            #metabNetwork_sif2.write('%s\tIdentified\t \t \t \n' % theMetab)

        if theMetab in metabFreq.keys():
            metabFreq[theMetab][0] += 1
        else:
            metabFreq[theMetab] = [1, 0, 0, metab_type]
        
    else:
        #metabGMT_file.write('%s\t%s\t' % (metabSetName, metabSetDescrip))
        for j in xrange(len(metabSet)):
            #if j == (len(metabSet)-1):
                #metabGMT_file.write('%s\n' % (metabSet[j].replace('"', '')))
            #else:
                #metabGMT_file.write('%s\t' % (metabSet[j].replace('"', '')))
            
            for k in xrange(j, len(metabSet)):
                if j==k:
                    continue

                metab1 = metabSet[j].replace('"', '').replace('X_', '')
                metab2 = metabSet[k].replace('"', '').replace('X_', '')
                #metabNetwork_sif.write('%s %s %s\n' % (metab1, model, metab2))

                
                try:
                    print(int(metab1))
                    metab1_ind = 2
                    metab1_type = "Unidentified"
                except ValueError:
                    if (metab1[0] == '_'):
                        pdb.set_trace()
                    metab1_ind = 1
                    metab1_type = "Identified"
                    

                try:
                    print(int(metab2))
                    metab2_ind = 2
                    metab2_type = "Unidentified"
                except ValueError:
                    if (metab2[0] == '_'):
                        pdb.set_trace()
                    metab2_ind = 1
                    metab2_type = "Identified"

                
                if metab1 in metabFreq.keys():
                    metabFreq[metab1][metab2_ind] += 1
                else:
                    if metab2_ind == 1:
                        metabFreq[metab1] = [0, 1, 0, metab1_type]
                    else:
                        metabFreq[metab1] = [0, 0, 1, metab1_type]

                if metab2 in metabFreq.keys():
                    metabFreq[metab2][metab1_ind] += 1
                else:
                    if metab1_ind == 1:
                        metabFreq[metab2] = [0, 1, 0, metab2_type]
                    else:
                        metabFreq[metab2] = [0, 0, 1, metab2_type]
                if (modelType == 'comparison'):
                    #if model in mainModels:
                    #    metabNetwork_sif2.write('%s\t%s\t%s\t%s\t%s\n' % (metab1, metab1_type, "Genetic", metab2, metab2_type))
                    #    metabNetwork_sif2.write('%s\t%s\t%s\t%s\t%s\n' % (metab2, metab2_type, "Genetic", metab1, metab1_type))
                            
                    #elif model in dietModels:
                    #    metabNetwork_sif2.write('%s\t%s\t%s\t%s\t%s\n' % (metab1, metab1_type, "GxD", metab2, metab2_type))
                    #    metabNetwork_sif2.write('%s\t%s\t%s\t%s\t%s\n' % (metab2, metab2_type, "GxD", metab1, metab1_type))

                    if (metab1 in all_degreeDist.keys()):
                            all_degreeDist[metab1] += 1
                    else:
                        all_degreeDist[metab1] = 1

                    if (metab2 in all_degreeDist.keys()):
                        all_degreeDist[metab2] += 1
                    else:
                        all_degreeDist[metab2] = 1
                else:
                    print('Hello!')
                    #metabNetwork_sif2.write('%s\t%s\t%s\t%s\t%s\n' % (metab1, metab1_type, model, metab2, metab2_type))
                    #metabNetwork_sif2.write('%s\t%s\t%s\t%s\t%s\n' % (metab2, metab2_type, model, metab1, metab1_type))
                

#metabNetwork_sif.close()
#metabNetwork_sif2.close()
#metabGMT_file.close()

metabFreq_file = open('Z:/Drosophila/Bayesian/general/results/metabolites/updated/plots/metabFreq_%s.csv'%modelType, 'w')
metabFreq_file.write('Metabolite,Itself,Other-ID,Other-NonID,MetabType\n')
for metab in metabFreq.keys():
    metabFreq_file.write('%s,%i,%i,%i,%s\n' % (metab, metabFreq[metab][0], metabFreq[metab][1], metabFreq[metab][2], metabFreq[metab][3]))
metabFreq_file.close()

networkDegreeDist_file = open('Z:/Drosophila/Bayesian/general/results/metabolites/updated/plots/networkDegreeDist_%s.csv'%modelType, 'w')
networkDegreeDist_file.write('Metabolite,DegreeDist\n')
for m in all_degreeDist.keys():
    networkDegreeDist_file.write('%s,%i\n' % (m, all_degreeDist[m]))
networkDegreeDist_file.close()




