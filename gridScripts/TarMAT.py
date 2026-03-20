import os
if __name__ == '__main__':
    #Tarring MAT opts folder
    cmd = "tar -cvzf /exp/minerva/app/users/alhart/opt.tar.gz -C /exp/minerva/app/users/alhart/MAT_AL9/opt/ ."
    os.system(cmd)

    #If we already have a tarred opts folder, remove it so I can copy over new one (/persistent/ doesn't like being overwritten)
    if (os.path.isfile("/pnfs/minerva/persistent/users/alhart/NuMuNukeIncl/TarredMATFramework/opt.tar.gz")):
        cmd = "rm /pnfs/minerva/persistent/users/alhart/NuMuNukeIncl/TarredMATFramework/opt.tar.gz"
        os.system(cmd)

    #Copy opts tar
    cmd = "mv /exp/minerva/app/users/alhart/opt.tar.gz /pnfs/minerva/persistent/users/alhart/NuMuNukeIncl/TarredMATFramework/opt.tar.gz"
    os.system(cmd)