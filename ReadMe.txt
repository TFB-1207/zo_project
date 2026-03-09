ReadMe

Structure 
The main folders are campaigns, manuscripts and presentations.
campaigns: data and meta-data is organized per campaign per site. All materials related to a specific site and campaign (e.g. Moehne winter campaign) are stored within a single directory. This is both fieldwork preparation files, such as permits, material lists etc., data produced from that campaign and scripts to process the data. The structure of these field campaigns directories are universal. 
data
raw data is stored under raw (only on /W drive?). Raw sequencing data are processed by project partners (Kevin) and are not stored locally. Partner-delivered data products are stored as preprocessed data in the folder structure (what we get from Kevin). All subsequent transformations performed within this project are stored as processed data.
The idea is that data in the campaigns folder is ultimately so that we have a neat overview from what is gathered in that campaign, to be later used when creating data for a chapter/publication. Therefore, manuscripts will have a seperate scripts and results folder. 
fieldwork
The field folder contains all meta-data aquired during the field work. the folder prep contains all information regarding the preperation of the fieldwork. 
scripts
contains all scripts to pre-process and process the raw data. 
manuscripts: This folder contains the scripts and data used in the actual chapter/publication. This can be aquired from multiple fieldwork campaigns. f.e. combining sequencing and sediment property data from both the Moehne winter and summer campaigns. 
It is structured per chapter per PhD condidate (Brian = BRC, Tom = TFB).
presentations:
contains figures, tables and ppx files used during poster presentations or talks. 
 

Naming conventions
Site abbreviations consist of the first two letters of the site name, with the first letter capitalized (e.g. Moehne -> Mo).
Campaigns are abbreviated by season: winter (wi), spring (sp), summer (su) or fall (fa).
Folder names follow the structure [site]_[season] and then within that folder you find fieldwork, data and scripts folders. An example at the end of this document. The actual files files will also follow the [site]_[season] structure (e.g. Mo_wi_grainsize_raw.csv)
field samples: [site abbreviation]_[season abbreviation]_[unique number] (e.g. Mo_wi_35)

Capslock usage
I prefer no caps unless it is for names or abbreviations (e.g. Moehne or TOC). This will save typos and frustrations.

Abbreviations
See xlsx file. 


Scripts
In SharePoint we only upload the scripts that produce the data also saved in SharePoint. F.e. the script which takes the raw grainsize data and produces the pre-processed data set. The idea is that, when setting the right working directory, a script in SharePoint should run fully. So all variables used are internally created in the script, from data which is available to us all. For this it is important to keep the site and fieldcampaign naming in the variables, otherwise everything will overwrite when we have multiple campaigns. e.g. df_plants_genus is universal, so this should be Mo_wi_df_pl_gs, this way the dataframe with plant genera from the Moehne winter samples does not get overwritten when creating a similar dataframe in another field campaign.  

---------------------------------------------------------------------------------------------------------------------------------------
folder structure example:
|RESIST project Documents
	-campaigns
		- Moehne
			- Mo_wi
 				- fieldwork
					- prep
						- samp_strat
							-Mo_wi_sampling_locations.tiff 
						- materials
							- Mo_wi_materials_list.docx
						- permits
							- Mo_wi_permits.docx
						- logistics
							- Mo_wi_logistics.docx
					- field
						- samples_records
							- Mo_wi_sample_list_field.xlsx 
						- field_photos
							-site_overview
								- Mo_wi_[date].jpg
							- sampling_points
								- Mo_wi_[sample].jpg
							- cores
								- Mo_wi_[core_name]_1.jpg
						- field_notes
							- Mo_wi_field_note_[date]_[# of note].txt
				- data
					- raw (raw data is only saved at W drive, so not in our sharepoint)
						- sed_raw
							- Mo_wi_sed_raw.csv
						- grainsize_raw
							- Mo_wi_grainsize_raw.csv
						- OSL_raw
						- Rock_Eval_raw
					- pre_processed
						- sediments
 							- Mo_wi_grainsize.csv
 							- Mo_wi_Rock_Eval.csv 
						- sequencing
							- Mo_wi_seq_[date].tsv
					- processed
						- sediments
							- Mo_wi_sed_overview.csv
						- sequencing
							- Mo_wi_gs_pl.tsv
					- results
						- figures
							- deamination plots
								- final
							- abundance 
							- diversity 
						- tables
						- statistics
				- scripts (only upload the final, working script which produces the according data, no work in progress 					   scripts). 
					- pre_processing
						- sediments
							- Mo_wi_grainsize_pre_pro.rmd
							- Mo_wi_rock_eval_pre_pro.rmd
						- sequence_data
							- Mo_wi_gs_pl.rmd
				
			- Moehne_su (this will contain an identical structure)
	- manuscripts 
		- BRC_chapter01
			- results
				- figures
				- tables
			- scripts
	- presentations
		- posters
		- talks

