# UniProt .dat file parser
Module for easily parsing data from UniProt's .dat files on their FTP server. Parses the .dat file for data regarding protein's name, isoforms, description, Accession ID, gene, and database(ie TrEMBL or SwissProt). First download the uniprot_dat_parser module and install on system, then follow these instructions:
  1. Download UniProt .dat file from the ftp server: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/ (This site is from the Taxonomic divisions file section on the ftp server)
  2. Download the uniprot_dat_parser.py module and ensure that the parser is in the same directory as the .dat uniprot file.
  3. Import the uniprot_dat_parser (import uniprot_dat_parser as udp)
  4. Parse the file using the parse function (uniprot_dict = udp.parse("uniprot_filename.dat"))
  5. Access information from the data dictionary created:
    a. To access isoform information, type (uniprot_dict[Accession_ID]["Isoform"]) \n
    b. To access Gene information, type (uniprot_dict[Accession_ID]["Gene"]) \n
    c. To access Description information, type (uniprot_dict[Accession_ID]["Description"]) \n
    d. To access Database information, type (uniprot_dict[Accession_ID]["DB"]) \n
    e. To access Sequence Variant information, type (uniprot_dict[Accession_ID]["SV"]) \n
    f. To access Protein Existence information (1-5), type (uniprot_dict[Accession_ID]["PE"]) \n
