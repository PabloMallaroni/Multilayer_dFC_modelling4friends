# Multilayer_dFC_modelling
Code serves the purpose to look at whole brain static modularity and dynamic FC integ/seg as ^suggested^ by 
https://www.nature.com/articles/s41591-022-01744-z

Feel try to try on any of your data for replication and check the code.
These scripts were written to account for potential roi numbers per subject.
Thus, works with cells of subject x session containing timeseries x roi data as input to step 1

The idea is to have consistent code for this kind of approach so please try on your stuff
If you have suggestions for better sliding windows please feel free to branch and debug

TBD: step 5 with plots and stats (depends on your design):
but you can treat step4  as a net x net ttest 
and step 2 q_norm result as a single array of sub x cond normalised q_modularity values
