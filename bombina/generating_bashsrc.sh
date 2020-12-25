#!/bin/bash -
#title           :generating_bashsrc.sh
#description     :This script will make a header for a bash script.
#author		       :bgw
#date            :20111101
#version         :0.4
#usage		       :bash mkscript.sh
#notes           :
#bash_version    :4.1.5(1)-release
#==============================================================================

today=$(date +%Y%m%d)
div=======================================

/usr/bin/clear

_select_title(){

    # Get the user input.
    printf "Enter a title: " ; read -r title

    # Remove the spaces from the title if necessary.
    title=${title// /_}

    # Convert uppercase to lowercase.
    title=${title,,}

    # Add .sh to the end of the title if it is not there already.
    [ "${title: -3}" != '.sh' ] && title=${title}.sh

    # Check to see if the file exists already.
    if [ -e $title ] ; then
        printf "\n%s\n%s\n\n" "The script \"$title\" already exists." \
        "Please select another title."
        _select_title
    fi

}

_select_title

printf "Enter a description: " ; read -r dscrpt
name="Ruby(Yiru) Sheng"
printf "Enter the version number: " ; read -r vnum

# Format the output and write it to a file.
printf "%-16s\n\
%-16s%-8s\n\
%-16s%-8s\n\
%-16s%-8s\n\
%-16s%-8s\n\
%-16s%-8s\n\
%-16s%-8s\n\
%-16s%-8s\n\
%-16s%-8s\n\
%s\n\n\n" '#!/bin/bash' '#title' ":$title" '#description' \
":${dscrpt}" '#author' ":$name" '#date' ":$today" '#version' \
":$vnum" '#usage' ":./$title" '#notes' ':' '#bash_version' \
":${BASH_VERSION}" \#$div${div} > $title

# Make the file executable.
chmod +x $title

/usr/bin/clear

_create_anotherfile(){
  printf "%s\n%s\n%s\n%s\n%s\n%s\n%s\n\n" "Which step do you want to process? (input with number please)" \
         "1 for trim / quality control." \
         "2 for assembly." \
         "3 for quantify." \
         "4 for generate expression matrix." \
         "5 for annotation." \
         "0 for other steps."
  read -r stepchosen
  echo $stepchosen
  #Select between Vim or Emacs.
      printf "%s\n%s\n%s\n" "If you want to start with another script with new header" "1 for Yes." "2 for No."
      read -r newfiledecision


      # Open the file with the cursor on the twelth line.
      case $newfiledecision in
          1) source generating_bashsrc.sh
              ;;
          2) exit
              ;;
          *) /usr/bin/clear
             printf "%s\n%s\n\n" "I did not understand your selection." \
                 "Press <Ctrl-c> to quit."
             _create_anotherfile
              ;;
      esac
}
_create_anotherfile
# _select_editor(){
#
#     # Select between Vim or Emacs.
#     printf "%s\n%s\n%s\n\n" "Select an editor." "1 for Vim." "2 for Emacs."
#     read -r editor
#
#     # Open the file with the cursor on the twelth line.
#     case $editor in
#         1) vim +12 $title
#             ;;
#         2) emacs +12 $title &
#             ;;
#         *) /usr/bin/clear
#            printf "%s\n%s\n\n" "I did not understand your selection." \
#                "Press <Ctrl-c> to quit."
#            _select_editor
#             ;;
#     esac
#
# }
#
# _select_editor
