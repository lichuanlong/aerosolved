#!/bin/bash

AEROSOLVEDROOT=../../
GITINFOFILE=aerosolModel/aerosolModelGitInfo.H
TEMPGITINFOFILE=tempGitInfoFile.txt

if [ -x "$(command -v git)" ] && [ -d $AEROSOLVEDROOT/.git ]; then

    # When generated, the git info file's timestamp is set to now minus 1
    # minute. If the file is older than now minus two minutes, regenerate.
    # Without this check, the git info file would be regenerated at each .C file
    # compilation step, which would make things slower.

    OLDFILE=$(find aerosolModel -maxdepth 1 -name aerosolModelGitInfo.H -mmin +2)

    if [ ! -f "$GITINFOFILE" ] || [ -f "$OLDFILE" ]; then

        # Remove if the file already exists

        if [ -f "$GITINFOFILE" ]; then

            rm -f $GITINFOFILE

        fi

        GITBRANCH=$(git --git-dir $AEROSOLVEDROOT/.git rev-parse --abbrev-ref HEAD)
        GITCOMMIT=$(git --git-dir $AEROSOLVEDROOT/.git rev-parse HEAD)

        echo "void Foam::aerosolModel::gitInfo()" >> $GITINFOFILE
        echo "{" >> $GITINFOFILE

        echo "    Info<< \"Git state at compilation:\" << nl << nl" \
            >> $GITINFOFILE

        echo "        << \"    Git branch = \" << \"${GITBRANCH}\" << nl" \
            >> $GITINFOFILE

        echo "        << \"    Git commit = \" << \"${GITCOMMIT}\" << nl;" \
            >> $GITINFOFILE

        # Write diffs temporarily to file

        git diff --stat=999 > $TEMPGITINFOFILE

        if [ "$(cat $TEMPGITINFOFILE)" != "" ]; then

            echo "    Info<< \"    Git diff =\" << nl;" >> $GITINFOFILE

            cat $TEMPGITINFOFILE | while read LINE; do

                echo "    Info<< \"        $LINE\" << nl;" >> $GITINFOFILE

            done

        else

            echo "    Info<< \"    No differences found w.r.t. this commit\" << nl;" \
                >> $GITINFOFILE

        fi

        echo "}" >> $GITINFOFILE

        # Remove temporary file

        if [ -f "$TEMPGITINFOFILE" ]; then

            rm $TEMPGITINFOFILE

        fi

        # Reset the timestamp of the git diff file, to avoid an infinite loop
        # between aerosolModel.C dependencies

        touch -r $GITINFOFILE -d '-1 minute' $GITINFOFILE
    fi

else

    OLDFILE=$(find aerosolModel -maxdepth 1 -name aerosolModelGitInfo.H -mmin +2)

    if [ ! -f "$GITINFOFILE" ] || [ -f "$OLDFILE" ]; then

        echo "void Foam::aerosolModel::gitInfo()" > $GITINFOFILE
        echo "{" >> $GITINFOFILE
        echo "    Info<< \"    Git is not installed or this is not a git repository\" << nl;" \
            >> $GITINFOFILE
        echo "}" >> $GITINFOFILE

        # Reset the timestamp of the git diff file, to avoid an infinite loop
        # between aerosolModel.C dependencies

        touch -r $GITINFOFILE -d '-1 minute' $GITINFOFILE

    fi
fi
