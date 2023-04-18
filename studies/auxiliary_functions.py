import os
import json
import re
import ntpath
import jellyfish


process = '''process mappingFastQ {
    /*
    Mapping the read to genome with STAR
    Inputs : the reads R1 and R2 and the genome dir
    Outputs : file bam aligned
    */

    label 'STAR'
    publishDir 'data/map', mode: 'copy'

    input:
        tuple val(sample), path(R1), path(R2)
        path GenomeDir

    output: 
        path "*.bam"

    script:
        """
        STAR --outSAMstrandField intronMotif \
            --outFilterMismatchNmax 4 \
            --outFilterMultimapNmax 10 \
            --genomeDir ${GenomeDir} \
            --readFilesIn ${R1} ${R2} \
            --runThreadN ${params.nb_threads_star} \
            --outSAMtype BAM SortedByCoordinate \
            --outStd BAM_SortedByCoordinate \
            > ${sample}.bam
        """
}'''

#Fonction qui retourne les éléments d'intersection entre 2 ensembles
def intersection(l1, l2):
    return list(set(l1) & set(l2))

#Fonction qui retourne les éléments d'union entre 2 ensembles
def union(l1, l2):
    return list(set(l1 + l2))

#Indice de Jaccard
def jaccard(t1, t2):
    num = len(intersection(t1, t2))
    if(num==0):
        return 0
    denum= len(union(t1, t2))
    return num/denum

def test_george():
    print("yo")

#Function that returns the terms of the functions
def get_terms_function(dico):
    tab=[]
    for d in dico:
        for o in d['operation']:
            tab.append(o['term'].lower())
    return list(set(tab))

#Function that returns the terms of the topics
def get_terms_topic(dico):
    tab=[]
    for d in dico:
        tab.append(d['term'].lower())
    return tab

#Functions that gets the bio.tools description, function terms and topics for a given bioinformatics tool 
def get_info(tool = 'star'):
    #--------------------------------------------------
    try:
        with open("bio_dot_tools_archive.json") as json_file:
            archive = json.load(json_file)
    except:
        archive = {}
    #--------------------------------------------------
    try:
        return archive[tool]
    except:
        try:
            command = f'curl -X GET "https://bio.tools/api/{tool.lower()}/?format=json" > temp.json'
            os.system(command)
            f = open("temp.json")
            os.system("rm temp.json")
            dico = json.load(f)
            res = {'description' : dico['description'].lower(), 'function' : get_terms_function(dico['function']), 'topic' : get_terms_topic(dico['topic'])}
            archive[tool] = res
            with open("bio_dot_tools_archive.json", "w") as outfile:
                json.dump(archive, outfile, indent = 4)
            return res
        except:
            #If nothing is found we save that too in the dictironnary
            archive[tool] = None
            with open("bio_dot_tools_archive.json", "w") as outfile:
                json.dump(archive, outfile, indent = 4)
            return None
    
#Function that returns the name of a given process
def get_name_process(process):
    result = re.search(r"process\s+([^\s\{]+)\s*{", process)
    return result.group(1)

#Function that finds the end of the definition of a long script in a process
def get_end_long_stript(process, start, patternMatch):
    for i in range(start, len(process)-2):
        if(process[i:i+3]==patternMatch[1:-1]): #patternMatch[1:-1] since with don't want the '()'
            return i
    raise Exception(f"Couldn't find the find of the strict in {process} with the patternMatch {patternMatch}")

#Function that extracts the multiple scripts from the given process
def get_multiple_scripts(process_input):
    process = process_input
    continue_searching = True
    scripts = []
    patternStart = [r'(""")', r"(''')"]
    while(continue_searching):
        #We are gonna be removing the patternStart as we go -> so if we can't find anymore we stop searching
        continue_searching = False
        for pattern in patternStart:
            for match in re.finditer(pattern, process):
                continue_searching = True
                patternMatch = pattern
                start = match.span(1)[1]
                end = get_end_long_stript(process, start, patternMatch)
                scripts.append(process[start:end].strip())
                process = process.replace(process[start-3:end+1+3], '')
                break #When we find one, we stop to start searching again from 0
    return scripts

#Function that removes the stub from the process -> we assume that it is the last thing defiened in the process 
def remove_stub(process):
    for match in re.finditer(r"stub *:\s", process):
        return process[:match.span()[0]].strip()+'\n}'
    return process

#Function that extracts the script(s) from a process
def get_script(process):
    #Remove stub
    process = remove_stub(process)
    #Case there are multiple scripts in the process
    if(process.count('"""')>2 or process.count("'''")>2):
        return get_multiple_scripts(process)
    else:
        patternStart = [r'(""")', r"(''')", r'\n+\s*"([^"]+)"\n', r"\n+\s*'([^']+)'\n", r'(script\s*:)\s*', r'(shell\s*:)\s*', r'(exec\s*:)\s*']
        for pattern in patternStart:
            found_match = False
            for match in re.finditer(pattern, process):
                patternMatch = pattern
                start = match.span(1)[1]
                group = match.group(1)
                found_match = True
                break
            if(found_match):
                break
        
        #Case it's a long script -> eg '"""' or "'''"
        if patternMatch == patternStart[0] or patternMatch == patternStart[1]:
            end = get_end_long_stript(process, start, patternMatch)
            return [process[start:end].strip()]
        
        #Case it's a single line script 
        elif patternMatch == patternStart[2] or patternMatch == patternStart[3]:
            return [group]
        
        #Case the script is defined using script, shell or exec -> we assume that the script end at the end of the definition of the process
        elif patternMatch == patternStart[4] or patternMatch == patternStart[5] or patternMatch == patternStart[6]:
            return [process[start:-1].strip()]
        
        else:
            return []

#Function that extarcts the language in a #!bin/.. commad
def extract_language(bash_command):
    if(len(bash_command.split())>1):
        return bash_command.split()[-1]
    return bash_command.split('/')[-1]

#Function that returns the languages in a script
def get_language(scripts):
    languages = []
    for s in scripts:
        found_match = False
        for match in re.finditer(r'#!([^\n]+)', s):
            found_match = True
            language = extract_language(match.group(1).strip())
            if language == 'sh' or language == 'ksh' or language == 'bashlog':
                language = 'bash'
            languages.append(language)
        if(not found_match):
            languages.append('bash')
    return languages


#Marines Code
#============================================

EXCEPTIONS = ["export", 'awk', 'sed', 'grep', "cmd", "module", "cat", "elif", "sort", "cd", "zcat",
              "rm", "for", "find", "java", "forgi", "sleep", "tabix", "zgrep", "wget", "mv", "mkdir", "echo",
              'FS', 'head', 'Rscript', "python", "jekyll","bgzip", "tr", "dot", "tRNA", "header", "fi", "then",
              "read", "do","else","cut", "wc", "tar", "gzip", "cool", "if", "turn", "git", "checkm", "cp", "make", "pour",
              "NR", "melt" , "read", "tail", "genes", "add", "bc", "scp", "scif", "uniq", "ln", "set", "zip", "time", "ls",
              "print","make", "pour", "source", "melt", "paste", "split", "layer", "touch" , "google-chrome", "query"
              "curl", "snps", "def"]

def path_leaf(path):
    """
    author : marinedjaffardjy 
    """
    # extract last word from path -- useful for extracting toolnames
    # input : string path name
    # output string last element of path
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

def parse_lines(compil):
    """
    author : marinedjaffardjy 
    """
    # returns a list of lines from the shell data
    # input : shell content (string)
    # output : list of command lines as a list of strings
    lines_list = []
    for whole_line in compil.split('\n'):
        # check if the line has at least one character
        if (re.search('[a-zA-Z]', whole_line) is not None):  # checks that the line contains a character
            # split the line into several lines if there is a |, a ; or a >
            if any(ext in whole_line for ext in ['', '|', ';', '>']):
                whole_line = whole_line.split('|')
                for subline in whole_line:
                    subline = subline.split('>')
                    for subsubline in subline:
                        subsubline = subsubline.split(';')
                        for el in subsubline:
                            lines_list.append(el)
    return lines_list

def get_toolname_from_line(line):
    """
    author : marinedjaffardjy 
    """
    # get command from a line
    # input : a line (string)
    # output : a command word (string) or a tuple of two words (string, string) or none if the toolname is not relevant
    while (len(line) > 1 and line[0] == ' '):
        line = line[1:]
    line = line.replace('cmd ( "', "")
    tool = line.split(' ')[0].replace('\t', '').replace('"', '').replace("'",
                                                                             '')  # just for syntax, take out the tabs
    if ("=" in tool):
        tool.split('=')[0]
    if ("'" in tool):
        tool = tool.replace("'", '')
    # if the tool is in a command with a path ex /home/marine\.outil
    if ('/' in tool):
        tool = path_leaf(tool)
    if ('bowtie-' in tool):
        tool = "bowtie"
    if ('htseq-count' in tool):
        tool = "htseqcount"
    # look for second word :
    if len(line.split(' ')) > 1:
        sec_word = line.split(' ')[1]  # weget the second word
        # if the second word starts with a letter and is bigger than one letter, we take it into account
        if (len(sec_word) > 1):
            if (re.search('[a-zA-Z]', sec_word[0]) is not None and len(sec_word) > 1):
                # if the second word doesn't have an extension, we take it into account
                if ('.' not in sec_word):
                    # print("second word " + sec_word)
                    return [tool, sec_word]
    return [tool]

def get_toolnames(strScript):
    """
    author : marinedjaffardjy 
    Add some little changes : 
    Insteed of working on all Wf -> just the string of the script
    """
    toolnames = []
    scripts_python = []
    scripts_R = []
    scripts_bash = []
      
    lines_list = parse_lines(strScript)
    for line in lines_list:
        # print("line "+line)
        toolname = get_toolname_from_line(line)
        # print(toolname)
        if (re.search('[a-zA-Z]', toolname[0]) is not None and len(
                    toolname[
                            0]) > 1):  # if the toolname is comprised of at least one letter and is of length sup to one
            if (re.search('[a-zA-Z]', toolname[0][
                0]) is not None):  # if the toolname doesn't start with a special character
                    # # on cherche les scripts python
                    # if toolname[0] == "python":
                    #     scripts_python.append(toolname)
                    # # on cherche els scripts R
                    # if toolname [0]== "Rscript":
                    #     scripts_R.append(toolname)
                    # # on cherche les scripts bash
                    # if toolname[0][-3:]== ".sh":
                    #     scripts_bash.append(toolname)

                    if (toolname[0] not in EXCEPTIONS ):
                        if ('(' not in toolname[0] and '{' not in toolname[0] and '#' not in toolname[0] and "=" not in toolname[0] and '\\' not in toolname[0] and '.' not in toolname[0]):
                            toolnames.append(toolname)
                            # print("line " + line)
    #print(toolnames)
    return toolnames

#============================================



#Clemences Code
#============================================
#ClÃ©mence SEBE

# InspirÃ© du travail de : Marine Djaffardjy : https://github.com/mdjaffardjy/Snakemake_workflow_analysis/tree/main/src 

EXCEPTIONS_CLEM = ["export", 'awk', 'sed', 'grep', "cmd", "module", "cat", "elif", "sort", "cd", "zcat",
              "rm", "for", "find", "java", "forgi", "sleep", "zgrep", "wget", "mv", "mkdir", "echo",           # "tabix" ??
              'FS', 'head', 'Rscript', "python", "jekyll","bgzip", "tr", "dot", "tRNA", "header", "fi", "then",
              "read", "do","else","cut", "wc", "tar", "gzip", "cool", "if", "turn", "git", "cp", "make", "pour",     #"checkm" 
              "NR", "melt" , "tail", "genes", "add", "bc", "scp", "scif", "uniq", "ln", "set", "zip", "time", "ls",
              "print","make", "pour", "source", "melt", "paste", "split", "layer", "touch" , "google-chrome", "query"
              "curl", "snps", "def", 

              'pandoc', 'pdflatex', 'done', 'perl', 'egrep', 'rev'
              
              'alias', 'at', 'apropos', 'aspell', 'autoexpect', 'bash', 'bunzip2', 'bzip2', 'chgrp', 'chmod', 'chown', 'cpio',
              'cron', 'crontab', 'chsh', 'cvs', 'date', 'dd', 'df', 'diff', 'dpkg', 'du', 'disown', 'eject', 'env', 'exit', 'export', 
              'except', 'fdisk', 'fg','file','finger','ftp','g++','gcc','gftp','groups','gvimdiff','gunzip','halt','hexdump','history',
              'hostname','id','ifconfig','info','init','iptables','iptraf','jobs','kill','killall','ldd','less','lsmod','lsof','look',
              'man','md5sum','mkfs','minicom','more','mount','netcat','netstat','nice','nm','objdump','openssl','passwd','ping','ps','pwd',
              'quota','rar','reboot','rename','route','rpm','rsync','screen','setenv','shutdown','ssh','su','sudo','tcpdump','tee','top',
              'traceroute','tac','ulimit','umount','uname','unset','unzip','unrar','uptime','useradd','userdel','usermod','vim','Vgcreate',
              'Vgdisplay','Vgs','Vgscan','vmstat','who','which','whoami','write','xargs','xev','xkill','xosview','yacc','yes','yum','yast','yast',
              
              'switch','case','break','string','null','trap','jar' ,'template','while', 'until', 'sh'   #error ? 
              ]


#Nettoyer le script : si une commande est Ã©crite sur plusieurs lignes ==> on la ramÃ¨ne Ã  une seule ligne
def clean(script):
    scriptTemp = script.replace('\\\n', '')
    return scriptTemp


#A partir d'un script : le couper en un ensemble de lignes (en fonction des retours Ã  la ligne et des autres caractÃ¨res spÃ©ciaux pouvant indiquer qu'il y ait un outil aprÃ¨s)
def getLines(script):
    lines = []
    #print(len(script.split("\n")))
    for ligne in script.split("\n"):

        if (re.search('[a-zA-Z]', ligne) is not None):
            subLines = re.split('\||;|&', ligne)
            for sub in subLines:
                tmp = sub.strip()
                if len(tmp) != 0:
                    lines.append(tmp)
    return lines


def getCandidats(ligne):
    tabCandidat = []
    sep = ligne.split()
    
    if len(sep) == 0:
        return []
    
    #travail sur le premier mot
    firstWord = sep[0]
    if (re.search('[a-zA-Z|\/]', firstWord[0]) is not None) and len(firstWord) > 1 : #rajouter que le mot doit faire plus d'une lettres ?

        #print(firstWord)
        if '/' in firstWord and firstWord[0:2] != '//':
            temp = path_leaf(firstWord)
        else: 
            temp = firstWord

        if (re.search('[\(|{|#|=|\\|\.|}|\)|\$|\/|:|!|\']', temp) is not None):
            
            return []

        if not temp in EXCEPTIONS_CLEM:
            tabCandidat.append(temp)
        
        if temp in EXCEPTIONS_CLEM : 
            return []
    else : 
        return []


    #travail sur le second mot
    if len(sep) > 1 :
        secondWord = sep[1]

        if secondWord[0] == '=' or secondWord[0:2] == '+=' or secondWord[0] == ':': #on est en presence d'une def de variable
            return []
        
        if (re.search('[a-zA-Z|\/]', secondWord[0]) is not None): 
            if (re.search('[\(|{|#|=|\\|\.|}|\)|\$|\/]', secondWord) is not None):
                return tabCandidat
            else : 
                tabCandidat.append(secondWord) 
        else : 
            return tabCandidat

    return tabCandidat

def findCandidat(script):
    candidats = []
    lines = getLines(script)

    for l in lines :
        c = getCandidats(l)
        #print(l)
        if c != []:
            candidats.append(c)
    #print(candidats)
    return candidats

#============================================



def flatten(l):
    return [item for sublist in l for item in sublist]


def get_tools(scripts, code = "marine"):
    tools = []
    for s in scripts:
        if(code=="marine"):
            tools += flatten(get_toolnames(s))
        elif(code=="clemence"):
            tools += flatten(findCandidat(s))
    for i in range(len(tools)):
        tools[i] = tools[i].lower()
    return list(set(tools))



#TODO -> FIND A BETTER WAY OF DOING THIS
def get_descriptions(tools):
    des = ""
    for t in tools:
        dict = get_info(t)
        try:
            des += dict['description']+'-'
        except:
            None
    return des

#TODO -> FIND A BETTER WAY OF DOING THIS
def get_functions(tools):
    fun = []
    for t in tools:
        dict = get_info(t)
        
        try:
            fun += dict['function']
        except:
            None
    return fun

#TODO -> FIND A BETTER WAY OF DOING THIS
def get_topics(tools):
    top = []
    for t in tools:
        dict = get_info(t)
        
        try:
            top += dict['topic']
        except:
            None
    return top


def normalised_levenshtein(txt1, txt2):
    if(txt1!='' or txt2!=''):
        score = jellyfish.levenshtein_distance(txt1,txt2)
        l = max(len(txt1),len(txt2))
        return (l-score)/l
    return 0


#TODO -> FIND A BETTER WAY OF DOING THIS
def get_biocontainers_tags(tools):
    top = []
    with open('./archives/biocontainers.json') as json_file:
        biocontainers = json.load(json_file)
    for t in tools:
        #Try finding tool
        try:
            top += biocontainers[t]["tags"]
        except:
            None
    return top

#TODO -> FIND A BETTER WAY OF DOING THIS
def get_biocontainers_tags_last_word(tools):
    top = []
    with open('./archives/biocontainers.json') as json_file:
        biocontainers = json.load(json_file)
    for t in tools:
        #Try finding tool
        try:
            top += biocontainers[t]["tags"]
        except:
            None

    temp = []
    for i in range(len(top)):
        last_word = top[i].split(":")[-1]
        if(last_word!=""):
            temp.append(last_word)
    return temp

#TODO -> FIND A BETTER WAY OF DOING THIS
def get_biocontainers_summary(tools):
    des = ""
    with open('./archives/biocontainers.json') as json_file:
        biocontainers = json.load(json_file)
    for t in tools:
        #Try finding tool
        try:
            des += "-".join(biocontainers[t]["summary"])+'-'
        except:
            None
    return des


def compare_processes(process1, process2, tool_extractor = "marine"):
    dico={}
    dico["name_similarity"] = normalised_levenshtein(get_name_process(process1), get_name_process(process2))
    dico["process_similarity"] = normalised_levenshtein(process1, process2)
    dico["process_similarity_no_white_spaces"] = normalised_levenshtein(process1.replace(' ', ''), process2.replace(' ', ''))
    script1 = get_script(process1)
    script2 = get_script(process2)
    #TODO -> find a better way of concatenating the script
    dico["script_similarity"] = normalised_levenshtein('-'.join(script1), '-'.join(script2))
    dico["script_similarity_no_white_spaces"] = normalised_levenshtein('-'.join(script1).replace(' ', ''), '-'.join(script2).replace(' ', ''))
    dico["language_similarity"] = jaccard(get_language(script1), get_language(script2))
    tools1 = get_tools(script1, code=tool_extractor)
    tools2 = get_tools(script2, code=tool_extractor)
    dico["tool_similarity"] = jaccard(tools1, tools2)
    #TODO : Need to find a better way to concatenate the description, functions, topics in the case of multiple tools
    dico["description_similarity"] = normalised_levenshtein(get_descriptions(tools1), get_descriptions(tools2))
    dico["description_similarity_no_white_spaces"] = normalised_levenshtein(get_descriptions(tools1).replace(' ', ''), get_descriptions(tools2).replace(' ', ''))
    dico['function_similarity'] = jaccard(get_functions(tools1), get_functions(tools2))
    dico['topics_similarity'] = jaccard(get_topics(tools1), get_topics(tools2))
    #BioContainers
    #TODO -> need to find better to associate
    dico['biocontainers_tags_similarity'] = jaccard(get_biocontainers_tags(tools1), get_biocontainers_tags(tools2))
    dico['biocontainers_tags_last_word_similarity'] = jaccard(get_biocontainers_tags_last_word(tools1), get_biocontainers_tags_last_word(tools2))
    dico["biocontainers_summary_similarity"] = normalised_levenshtein(get_biocontainers_summary(tools1), get_biocontainers_summary(tools2))
    dico["biocontainers_summary_similarity_no_white_spaces"] = normalised_levenshtein(get_biocontainers_summary(tools1).replace(' ', ''), get_biocontainers_summary(tools2).replace(' ', ''))
    return dico
