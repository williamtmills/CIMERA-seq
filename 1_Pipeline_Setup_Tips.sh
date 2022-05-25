# Download VirtualBox (https://www.virtualbox.org/)
#	OS X hosts (https://download.virtualbox.org/virtualbox/6.1.22/VirtualBox-6.1.22-144080-OSX.dmg)

# Download Ubuntu (https://www.ubuntu.com/)
#	Download > Ubuntu Desktop > 20.04 LTS (https://ubuntu.com/download/desktop/thank-you?version=20.04.2.0&architecture=amd64)

# Change privacy settings on computer to run VirtualBox
#	Apple Logo > System Preferences > Security & Privacy
#	"System software from developer "Oracle America Inc." was blocked from loading. > Allow


# 1. Open VirtualBox
# 2. Select 'New'
# 3. Fill out the 'Name and operating system' form

#	Name: Ubuntu 2021
#	Machine Folder: /Users/williamtmillsiv/VirtualBox VMs
#	Type: Linux
#	Version: Ubuntu (64-bit)

# 4. Fill out the 'Memory size' form

#	8 GB (8192 MB) (This appears to be what is necessary to index mouse genome with HISAT2)

# 5. Fill out the 'Hard disk' form

#	Create a virtual hard disk now

# 6. Fill out 'Hard disk file type' form

#	VDI (VirtualBox Disk Image)

# 7. Fill out 'Storage on physical hard disk' form

#	Dynamically allocated

# 8. Fill out the 'File location and size' form

#	/Users/williamtmillsiv/VirtualBox VMs/Ubuntu 2021/Ubuntu 2021.vdi

#	20.00 GB

# 9. Click 'Create'
# 10. Click the 'Settings (Gear)' icon

#	System > Processor > 1 CPU (never more than half)

#	Display > Video Memory > 128 MB (maximum)

#	Display > Acceleration > Enable 3D Acceleration

#	Storage > Controler: IDE > Empty > Optical Drive > Click the 'Disk' icon > Choose/Create a Virtual Optical Disk > Add > Select 'ubuntu-20.04.2.0-desktop-amd64.iso' from Downloads folder

# 11. Create a new folder on your computer called 'Ubuntu_Sharing_2022'*

#	Shared Folders > Click the 'New Folder' icon

#	Folder Path: /Users/williamtmilsiv
#	Folder name: Ubuntu_Sharing_2022
#	Check 'Auto-mount

#	Click 'OK

# 12. Start Ubuntu 2021

# 13. Click 'Install Ubuntu'

#	Keyboard layout
#	Choose your keyboard layout: English (US) > English (US) > Continue

#	Updates and other software
#	What apps would you like to install to start with? > Minimal Installation (Web browser and basic utilities.)
#	Other options > Download updates while installing Ubuntu
#	Other options > Install third-party software for graphics and Wi-Fi hardware and additional media formats
#	Continue

#	Installation type
#	This computer currently has not detected operating system. What would you like to do? > Erase disk and install Ubuntu > Install Now > Continue

#	Where are you? > New York > Continue

#	Who are you?
#	Your name: William Mills
#	Your computer's name (autofill): william-VirtualBox
#	Pick a username (autofill): william
#	Choose a password: 
#	Confirm your password: 
#	Require my password to log in
#	Continue

# 14. Mount shared folder on virtual computer

#	Configure your Ubuntu workspace by running these commands in terminal

#	sudo apt-get update
#	sudo apt-get upgrade
#	sudo apt-get install build-essential gcc make perl dkms

#	Devices > Insert Guest Additions CD image... > Run > Restart Virtual Computer

#	If there is an error, try

#	File > Virtual Media Manager > Optical Disks > Select 'VBoxGuestAdditions.iso' > Click "Release" button > Click "Remove" button

#	sudo adduser william vboxsf

#	Close and re-open virtual computer

# 15. Download the appropriate version of miniconda

#	Verion of python installed in Ubuntu can be found by running:	python3 --version
#	Different versions of Miniconda for Linux can be found here: https://docs.conda.io/en/latest/miniconda.html#linux-installers

#	Run these commands in terminal

#	wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.11.0-Linux-x86_64.sh
#	bash Miniconda3-py38_4.11.0-Linux-x86_64.sh

#	Set path to Miniconda to simplify command-line usage
#	May have to edit based on where Miniconda is installed
#	If you used the installation above from your home directory, you should just have to change 'william' to your username

#	Run this command in terminal

#	PATH=$PATH:/home/william/miniconda3/condabin/

# 16. Create a conda environment where requisite packages are installed

#	Run this command in terminal
# 	Ensure the 'environment.yml file is in your shared folder

#	conda env create -f /media/sf_Ubuntu_Sharing_2022/environment.yml -n PACeR

