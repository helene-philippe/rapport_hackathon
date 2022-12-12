Pour executer la pipeline d'analyse : 

1 - activer l'environnement 'snakemake' avec la commande suivante : 

conda activate snakemake 

2 - Lancer le snakefile avec la commande suivante (remplacer <n> par le nombre de coeurs souhaités) :

snakemake -s snakefile --use-singularity --cores <n>


