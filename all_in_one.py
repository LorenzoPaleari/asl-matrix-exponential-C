import subprocess
import os

def roofline(version):

    values = [32, 128, 512, 1024]

    for i in values:
        for j in [0, 1, 2, 5]:
            file = open("src/benchmark/main.cpp", "r")
            lines = file.readlines()
            file.close()

            for k in range(len(lines)):
                if 'for (int t = ' in lines[k]:
                    lines[k] = 'for (int t = '+str(j)+'; t < '+ str(j+1)+'; t++){\n'

            file = open("src/benchmark/main.cpp", "w")
            file.writelines(lines)
            file.close()

            command = "make roofline N=" + str(i) + " V=" + version + " FLOP=N ROOFLIN=E"
            if "no_blas" in version or "mmm" in version or "solve" in version:
                command = "make roofline N=" + str(i) + " V=" + version + " FLOP=N ROOFLIN=E NO_BLA=S"

            proc = subprocess.run(command.split(), capture_output=True, text=True)
            print(proc.stdout)

            command = "mkdir data_collection/" + version + "/roofline"
            proc = subprocess.Popen(command.split(), stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
            output, _ = proc.communicate()

            command = "mv MyResults/ data_collection/"+ version +"/roofline/" + str(i) + "_" + str(6 - j)
            proc = subprocess.Popen(command.split(), stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
            output, _ = proc.communicate()

    file = open("src/benchmark/main.cpp", "r")
    lines = file.readlines()
    file.close()

    for k in range(len(lines)):
        if 'for (int t = ' in lines[k]:
            lines[k] = 'for (int t = 0; t < test_num; t++){\n'

    file = open("src/benchmark/main.cpp", "w")
    file.writelines(lines)
    file.close()


def adjust_register(version):
    file = open("src/benchmark/register.cpp", "r")
    lines = file.readlines()
    file.close()

    if "no_blas" in version or "mmm" in version or "solve" in version:
        for i in range(len(lines)):
            if 'add_function(&mexp,' in lines[i]:
                lines[i] = '  //add_function(&mexp, "mexp");\n'
            if 'add_function(&mexp_model' in lines[i]:
                lines[i] = '  //add_function(&mexp_model, "mexp_model");\n'
            if 'add_function_no' in lines[i]:
                lines[i] = '  add_function_no_blas(&mexp_no_blas, "mexp_no_blas");\n'
    elif "base" in version:
        for i in range(len(lines)):
            if 'add_function(&mexp,' in lines[i]:
                lines[i] = '  //add_function(&mexp, "mexp");\n'
            if 'add_function(&mexp_model' in lines[i]:
                lines[i] = '  add_function(&mexp_model, "mexp_model");\n'
            if 'add_function_no' in lines[i]:
                lines[i] = '  //add_function_no_blas(&mexp_no_blas, "mexp_no_blas");\n'
    else:
        for i in range(len(lines)):
            if 'add_function(&mexp,' in lines[i]:
                lines[i] = '  add_function(&mexp, "mexp");\n'
            if 'add_function(&mexp_model' in lines[i]:
                lines[i] = '  //add_function(&mexp_model, "mexp_model");\n'
            if 'add_function_no' in lines[i]:
                lines[i] = '  //add_function_no_blas(&mexp_no_blas, "mexp_no_blas");\n'
            
    file = open("src/benchmark/register.cpp", "w")
    file.writelines(lines)
    file.close()


def flop_run(version):
    command = "./test_runner_flop.sh"
    proc = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
    if "no_blas" in version or "mmm" in version or "solve" in version:
        output, _ = proc.communicate(input=version + " y")
    else:
        output, _ = proc.communicate(input=version + " n")

def clean_flops(version):
    file = open("data_collection/"+ version + "/flops.txt", "r")
    lines = file.readlines()
    file.close()

    values = [4, 8, 16, 32, 64, 128, 256, 512, 1024]
    indexes = []
    for i in values:
        for j in range(len(lines)):
            if "./bin " + str(i) in lines[j]:
                indexes.append(j)
                break
    
    flops = []
    file = open("data_collection/"+ version + "/flops_clean.txt", "w")
    bin = 0
    for index in indexes:
        return_index = []
        for i in range(1,7):
            for j in range(index,len(lines)):
                if "Performing test " + str(i) + "/6" in lines[j]:
                    for l in range(j, len(lines)):
                        if "> Checks " in lines[l]:
                            return_index.append(l)
                            break
                    break
        
        temp_flops = []
        k = 0
        for test in return_index:
            temp_in = [[0], [0], [0], [0]]
            for i in range(test + 1, len(lines)):
                if lines[i][0:2] != "Fl" and lines[i][0:2] != "NO" and lines[i][0:2] != "MA" and lines[i][0:2] != "UT" and lines[i][0:4] != "BASE":
                    break

                if "Flops:" in lines[i]:
                    temp_in[0][0] += float(lines[i].split()[1])
                elif "Flops N:" in lines[i]:
                    temp_in[1][0] += float(lines[i].split()[2])
                elif "Flops N^2:" in lines[i]:
                    temp_in[2][0] += float(lines[i].split()[2])
                elif "Flops N^3:" in lines[i]:
                    temp_in[3][0] += float(lines[i].split()[2])
            temp_flops.append(temp_in)
            k += 1
        bin += 1
        flops.append(temp_flops)

    for i in range(6):
        file.write("Test " + str(6 - i) + "\n")
        for j in range (0, 9):
            file.write("  N = " + str(values[j]) + "     ---     "+ str(flops[j][i][3][0]) + "N^3 " + str(flops[j][i][2][0])+ "N^2 "+ str(flops[j][i][1][0]) + "N " + str(flops[j][i][0][0])+"\n")
    file.write("\n\n\n")


def main():
    versions_to_test = ["solve_system_base",
                        "solve_system_base_opt",
                        "solve_system_base_vect",
                        "solve_system_lu",
                        "solve_system_lu_block",
                        "solve_system_lu_opt",
                        "solve_system_lu_vect",
                        "blas_mmm",
                        "base",
                        "base_no_blas",
                        "lu",
                        "lu_no_blas",
                        "opt1_basic",
                        "opt1_basic_no_blas",
                        "opt2_basic",
                        "opt2_basic_no_blas",
                        "opt3_basic",
                        "opt3_basic_no_blas",
                        "opt4_inline_basic",
                        "opt4_inline_basic_no_blas",
                        "opt4_inline_optimized",
                        "opt5_blocking_basic",
                        "opt5_blocking_basic_no_blas",
                        "opt6_vect_basic",
                        "opt6_vect_basic_no_blas",
                        "opt6_vect_lu",
                        "opt6_vect_lu_no_blas",
                        "opt6_vect_lu_not_lu",
                        "opt6_vect_lu_not_lu_no_blas",
                        ]
    
    for version in versions_to_test:
        command = "make cleanall"
        proc = subprocess.run(command.split(), capture_output=True, text=True)
        command = "mkdir data_collection/" + version
        proc = subprocess.run(command.split(), capture_output=True, text=True)

        adjust_register(version)

        roofline(version)
        flop_run(version)
        clean_flops(version)

        command = "./test_runner_data_collection.sh"
        proc = subprocess.Popen(command.split(), stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
        if "no_blas" in version or "mmm" in version or "solve" in version:
            output, _ = proc.communicate(input=version + " y")
        else:
            output, _ = proc.communicate(input=version + " n")

        command = "cp -r "+os.getcwd()+"/data/ data_collection/" + version + "/"
        proc = subprocess.Popen(command.split(), stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)


if __name__ == '__main__':
    main()
