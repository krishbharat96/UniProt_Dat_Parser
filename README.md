# UniProt dat file parser 
[![DOI](https://zenodo.org/badge/116890519.svg)](https://zenodo.org/badge/latestdoi/116890519)

Module for easily parsing data from UniProt's .dat files on their FTP server primarily for large UniProt database dumps. Parses the .dat file for data regarding protein's name, isoforms, description, Accession ID, gene, and database(ie TrEMBL or SwissProt). First download the uniprot_dat_parser module and install on system, then follow these instructions:
  1. Download UniProt .dat file from the ftp server: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/ (This site is from the Taxonomic divisions file section on the ftp server)
  2. Download the uniprot_dat_parser.py module and ensure that the parser is in the same directory as the .dat uniprot file:
  ```
  wget https://raw.githubusercontent.com/krishbharat96/UniProt_Dat_Parser/master/uniprot_dat_parser.py
  ```
  3. Import the uniprot_dat_parser 
  ```
  import uniprot_dat_parser as udp
  ```
  4. Parse the file using the parse function 
  ```
  uniprot_dict = udp.parse("uniprot_filename.dat")
  ```
  5. Access information from the data dictionary created:
  
        a. To access isoform information
        ```
        uniprot_dict[Accession_ID]["Isoform"]
        ```
        b. To access Gene information 
        ```
        uniprot_dict[Accession_ID]["Gene"]
        ```
        c. To access Description information 
        ```
        uniprot_dict[Accession_ID]["Description"]
        ```
        d. To access Database information (i.e. from SwissProt or TrEMBL)
        ```
        uniprot_dict[Accession_ID]["Database"]
        ```
        e. To access Sequence Variant information 
        ```
        uniprot_dict[Accession_ID]["SV"]
        ```
        f. To access Protein Existence information (1-5)
        ```
        uniprot_dict[Accession_ID]["PE"]
        ```
        g. To access Organism information
        ```
        uniprot_dict[Accession_ID]["Organism"]
        ```
        h. To access protein's ChEMBL ID
        ```
        uniprot_dict[Accession_ID]["Chembl_id"]
        ```
