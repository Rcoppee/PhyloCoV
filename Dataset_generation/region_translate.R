# Region translate
# antoine.bridier-nahmias@inserm.fr
# 2022/04/07

# These functions solve the problem of badly written region names in all files
# PAS MIS A JOUR !!! USE SHORT FUNCTION.
region_translate_long <- function(regions) {
  regions_all <- 
    list(
      "Grand-Est" = 'Grand-Est',
      "Grand Est" = 'Grand-Est',
      "Haut-Rhin" = 'Grand-Est',
      "Lorraine" = 'Grand-Est',
      "Champagne-Ardenne" = 'Grand-Est',
      "Nouvelle-Aquitaine" = 'Nouvelle-Aquitaine',
      "Nouvelle Aquitaine" = 'Nouvelle-Aquitaine',
      "Aquitaine" = 'Nouvelle-Aquitaine',
      "Auvergne-Rhone-Alpes" = 'Auvergne-Rhone-Alpes',
      "Auvergne et Rhône-Alpes" = 'Auvergne-Rhone-Alpes',
      "Haute-Savoie" = 'Auvergne-Rhone-Alpes',
      "Bourgogne-Franche-Comte" = 'Bourgogne-Franche-Comte',
      "Bourgogne et Franche-Comté" = 'Bourgogne-Franche-Comte',
      "Burgundy" = 'Bourgogne-Franche-Comte',
      "Bretagne" = 'Bretagne',
      "Centre-Val-de-Loire" = 'Centre-Val-de-Loire',
      "Centre-Val-De-Loire" = 'Centre-Val-de-Loire',
      "Centre-Val-De-Loire" = 'Centre-Val-de-Loire',
      "Centre-Val de Loire" = 'Centre-Val-de-Loire',
      "Corse" = 'Corse',
      "Hauts-de-France" = 'Hauts-de-France',
      "Nord-Pas-de-Calais" = 'Hauts-de-France',
      "Picardie" = 'Hauts-de-France',
      "Île-de-France" = 'Ile-de-France',
      "Ile-de-France" = 'Ile-de-France',
      "Ile-De-France" = 'Ile-de-France',
      "Occitanie" = 'Occitanie',
      "Normandie" = 'Normandie',
      "Normandy" = 'Normandie',
      "Basse-Normandie" = 'Normandie',
      "Pays-de-la-Loire" = 'Pays-de-la-Loire',
      "Pays de la Loire" = 'Pays-de-la-Loire',
      "Provence-Alpes-Cote-d'Azur" = 'Provence-Alpes-Cote-d-Azur',
      "Provence-Alpes-Côte-d'Azur" = 'Provence-Alpes-Cote-d-Azur',
      "Provence-Alpes-Côte d'Azur" = 'Provence-Alpes-Cote-d-Azur',
      "Provence-Alpes-Cote-d’Azur" = 'Provence-Alpes-Cote-d-Azur',
      "Marseille" = 'Provence-Alpes-Cote-d-Azur'
    )
  out_regions <- vector('list', length(regions))
  for (i in 1:length(regions)) {
    if (is.na(regions[i])) { 
      out_regions[[i]] <- "Unknown"
      next()      
    }
    out_regions[[i]] <- regions_all[[regions[i]]]
  }
  unlist(out_regions)
}

# This function solves the problem of badly written region names in all files
region_translate_short <- function(regions) {
  regions_all <- 
    list(
      "Grand-Est" = 'GES',
      "Grand Est" = 'GES',
      "Haut-Rhin" = 'GES',
      "Lorraine" = 'GES',
      "Champagne-Ardenne" = 'GES',
      "Nouvelle-Aquitaine" = 'NAQ',
      "Nouvelle Aquitaine" = 'NAQ',
      "Aquitaine" = 'NAQ',
      "Auvergne-Rhone-Alpes" = 'ARA',
      "Auvergne et Rhône-Alpes" = 'ARA',
      "Haute-Savoie" = 'ARA',
      "Bourgogne-Franche-Comte" = 'BFC',
      "Bourgogne-Franche-Comté" = 'BFC',
      "Bourgogne et Franche-Comté" = 'BFC',
      "Burgundy" = 'BFC',
      "Bretagne" = 'BRE',
      "Centre-Val-de-Loire" = 'CVL',
      "Centre-Val-De-Loire" = 'CVL',
      "Centre-Val-De-Loire" = 'CVL',
      "Centre-Val de Loire" = 'CVL',
      "Centre-val-de-Loire" = 'CVL',
      "Corse" = '20R',
      "Hauts-de-France" = 'HDF',
      "Nord-Pas-de-Calais" = 'HDF',
      "Picardie" = 'HDF',
      "Île-de-France" = 'IDF',
      "Ile-de-France" = 'IDF',
      "Ile-De-France" = 'IDF',
      "Occitanie" = 'OCC',
      "Normandie" = 'NOR',
      "Normandy" = 'NOR',
      "Basse-Normandie" = 'NOR',
      "Pays-de-la-Loire" = 'PDL',
      "Pays de la Loire" = 'PDL',
      "Provence-Alpes-Cote-d'Azur" = 'PAC',
      "Provence-Alpes-Côte-d'Azur" = 'PAC',
      "Provence-Alpes-Côte d'Azur" = 'PAC',
      "Provence-Alpes-Cote-d’Azur" = 'PAC',
      "Provence-Alpes-Cote-d-Azur" = 'PAC',
      "Marseille" = 'PAC'
    )
  out_regions <- vector('list', length(regions))
  for (i in 1:length(regions)) {
    if (is.na(regions[i])) { 
      out_regions[[i]] <- "Unknown"
      next()      
    }
    out_regions[[i]] <- regions_all[[regions[i]]]
  }
  unlist(out_regions)
}