#include <iostream>
#include <vector>
#include <numeric>
#include <string>
#include <fstream>
#include <sstream>
#include <queue>

using namespace std;
void lectura(vector<int> doctors_schedule_ref[5], int &M_maquinas_o, int &N_pacientes_o, int &D_doctores_o, int &P_r_o, int &P_p_o, int &P_u_o);
void paciente_queue(queue<int> &Q_paciente, vector<int> &L_pacientes, int day,  int P_r, int P_p, int P_u);
void DList(vector<int> doctor_schedule[5], int Lun, int Ma, int Mie, int Jue, int Vie);
void doctor_queue(queue<int> &Q_doctor, vector<int> &doctors_schedule_ref, int bloque);
void write(int schedule_pac_doctor[][20*16], int schedule_pac_machine[][20*16], int N_pacientes);
void PList(vector<int> &L_pacientes, int P_u, int P_p, int P_r);
void maquina_queue(queue<int> &Q_Maq, int M_maquinas);
void init_zeros(int matrix[][20*16], int N);
void clear(queue<int> &q );
void Prettywrite(int schedule_pac_doctor[][20*16], int schedule_pac_machine[][20*16], int N_pacientes);

int main() {
    int M_maquinas=0; int N_pacientes=0; int D_doctores=0; 
    int P_r=0; int P_p=0; int P_u=0;
    vector<int> doctors_schedule_ref[5];
    lectura(doctors_schedule_ref, M_maquinas, N_pacientes, D_doctores, P_r, P_p, P_u);
    cout << M_maquinas <<' '<<D_doctores <<' '<<N_pacientes << endl;
    cout << P_r <<' '<<P_p <<' '<<P_u << endl;
    // cout << doctors_schedule_ref[0][0] << endl;

    int days = 20; // 4semanas X 5 días
    int bloques = 16; //30 min each one
    
    int schedule_pac_doctor[N_pacientes][20*16];
    int schedule_pac_machine[N_pacientes][20*16];
    init_zeros(schedule_pac_doctor, N_pacientes);
    init_zeros(schedule_pac_machine, N_pacientes);
    
    vector<int> L_pacientes;

    queue<int> Q_paciente;
    queue<int> Q_Maq;
    queue<int> Q_doctor;
    int paciente_asigned; int doctor_asigned; int machine_asigned;

    PList(L_pacientes, P_u, P_p, P_r); // [N_sesions_p1, ...]

    for(size_t day = 0; day < days; day++){
        if(accumulate(L_pacientes.begin(), L_pacientes.end(),0)==0){break;}
        paciente_queue(Q_paciente, L_pacientes, day, P_r,P_p,P_u);
        for (size_t bloque = 0; bloque < bloques; bloque++) {
            // cout << " bloque "<< bloque << endl;
            maquina_queue(Q_Maq, M_maquinas);
            // cout << " en maquina " <<  Q_Maq.front();
            // cout << " con doctor " << Q_doctor.front() << endl;
            doctor_queue(Q_doctor, doctors_schedule_ref[day - int(day/5)*5], bloque);

            while(!Q_doctor.empty() && !Q_Maq.empty() && !Q_paciente.empty()){
                // cout << "dia "<< day+1 << " paciente "<<  Q_paciente.front()<< " en maquina "<<  Q_Maq.front()<< " con doctor " << Q_doctor.front() << endl;
                paciente_asigned = Q_paciente.front();
                doctor_asigned = Q_doctor.front();
                machine_asigned = Q_Maq.front();
                L_pacientes[paciente_asigned] -=1;
                Q_paciente.pop(); 
                Q_doctor.pop(); 
                Q_Maq.pop(); 

                schedule_pac_doctor[paciente_asigned][day*bloques + bloque] = doctor_asigned;
                schedule_pac_machine[paciente_asigned][day*bloques + bloque] = machine_asigned;
        
            if (accumulate(L_pacientes.begin(), L_pacientes.end(),0)==0) {
                cout << "break" << endl;
                break;
            }
            }
        }
    }
    
    Prettywrite(schedule_pac_doctor, schedule_pac_machine, N_pacientes);
    return 0;
}

void init_zeros(int matrix[][20*16], int N) {
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < 20*16; j++)
        {
            matrix[i][j]=0;
        }  
    }
}

void doctor_queue(queue<int> &Q_doctor, vector<int> &doctors_schedule_ref, int bloque) {
    /*
        doctors_schedule_ref: list of doctors availability -> [2,0,1,0,3] each field position represents one doctor
        bloque: number -> in [0,15]
    */
    clear(Q_doctor);
    for (size_t doctor = 0; doctor < doctors_schedule_ref.size(); doctor++) {
        if(doctors_schedule_ref[doctor] == 0) {
            continue;
        } else if(doctors_schedule_ref[doctor] == 1 && bloque < 8) { // 0: 9:00, 1: 9:30, 2: 10:00, 3: 10:30, 4: 11:00, 5: 11:30, 6: 12:00, 7: 12:30
            Q_doctor.push(doctor+1);
        } else if(doctors_schedule_ref[doctor] == 2 && bloque > 7) { // 8: 14:00, 9: 14:30, 10: 15:00, 11: 15:30, 12: 16:00, 13: 16:30, 14: 17:00, 15: 17:30
            Q_doctor.push(doctor+1);
        } else if(doctors_schedule_ref[doctor] == 3) { // ALL
            Q_doctor.push(doctor+1);
        }
   }

}

void clear(queue<int> &q ) {
   queue<int> empty;
   swap( q, empty );
}

void maquina_queue(queue<int> &Q_Maq, int M_maquinas) {
    clear(Q_Maq);
    for (size_t Maq = 0; Maq < M_maquinas; Maq++) {
        Q_Maq.push(Maq+1);
    }
}

void paciente_queue(queue<int> &Q_paciente, vector<int> &L_pacientes, int day,  int P_r, int P_p, int P_u) {
    clear(Q_paciente);
    for (size_t paciente = 0; paciente < L_pacientes.size(); paciente++){
        if(L_pacientes[paciente] != 0) {
           if(L_pacientes[paciente] !=2 && L_pacientes[paciente] !=4 && L_pacientes[paciente] !=30) {
                Q_paciente.push(paciente);
                // atendidos.append(paciente)
           } else{
                if(paciente < P_u){ // evalua prioridad urgente, paleativo, radical (paciente es un indice y la lisa que lo contiene esta ordenada segun tipos de pacientes)
                    // waiting_days = 2 - day
                    Q_paciente.push(paciente);

                } else if (paciente < P_u + P_p && 2 <= day+1){
                    // waiting_days = 14 - day
                    Q_paciente.push(paciente);

                } else if(paciente < P_u + P_p + P_r && 14 <= day+1 && (day+1)%5 != 0){
                    // waiting_days = 28 - day
                    Q_paciente.push(paciente);
                }
            }
        }   
    }
}


void PList(vector<int> &L_pacientes, int P_u, int P_p, int P_r) {
    for (size_t i = 0; i < P_u+P_p+P_r; i++) {
        if(i < P_u){
            L_pacientes.push_back(2);
        } else if (i < P_u + P_p) {
            L_pacientes.push_back(4);
        } else if (i < P_u + P_p + P_r) {
            L_pacientes.push_back(30);
        }
    }
}

void write(int schedule_pac_doctor[][20*16], int schedule_pac_machine[][20*16], int N_pacientes) {
    ofstream myfile;
    myfile.open("output");
    myfile << "           |                                          Day                                                         |\n";
    myfile << "           | B1 || B2 || B3 || B4 || B5 || B6 || B7 || B8 || B9 || B10 || B11 || B12 || B13 || B14 || B15 || B16 ||\n";
    for (size_t paciente_i = 0; paciente_i <  N_pacientes; paciente_i++) {
        myfile << " paciente" << paciente_i << ": " ;
        for (size_t doctor = 0; doctor < 20*16; doctor++)
        {
            myfile << " " << schedule_pac_doctor[paciente_i][doctor] << " ||";
        }
    myfile << "\n";
        
    }
    myfile.close();
}

void Prettywrite(int schedule_pac_doctor[][20*16], int schedule_pac_machine[][20*16], int N_pacientes) {
    ofstream myfile;
    myfile.open("output3");
    for (size_t day = 0; day < 20; day++)
    {    
        myfile << "           |                                          Day"<< day+1 <<"                                                       |\n";
        myfile << "           | B1 || B2 || B3 || B4 || B5 || B6 || B7 || B8 || B9 || B10 || B11 || B12 || B13 || B14 || B15 || B16 ||\n";
        for (size_t paciente_i = 0; paciente_i <  N_pacientes; paciente_i++) {
            myfile << " paciente" << paciente_i << ": " ;
            for (size_t doctor = day*16; doctor < (day+1)*16; doctor++)
            {
                myfile << " " << schedule_pac_doctor[paciente_i][doctor] << "  ||";
            }
        myfile << "\n";
            
        }
    }
    myfile.close();
}
void lectura(vector<int> doctors_schedule_ref[5], int &M_maquinas_o, int &N_pacientes_o, int &D_doctores_o, int &P_r_o, int &P_p_o, int &P_u_o) {
    ifstream archivo;
    int M_maquinas=0; int D_doctores=0; ; int N_pacientes=0;
    int Lun;int Ma;int Mie;int Jue;int Vie;
    int P_r=0; int P_p=0; int P_u=0;
    archivo.open("Caso4", ios::in);
    if(archivo.fail()){
        cout<<"NO se pudo"<<endl;
        exit(1);
    }

    int c = 0;
    for( string line; getline( archivo, line ); ) {
        istringstream iss(line);
        if (c==0) {
            iss >> M_maquinas; iss >> D_doctores; iss >> N_pacientes;
            // cout << M_maquinas <<' '<<D_doctores <<' '<<N_pacientes << endl;

        } else if (c<=D_doctores) {
            iss >> Lun; iss >> Ma; iss >> Mie; iss >> Jue; iss >> Vie;
            DList(doctors_schedule_ref, Lun, Ma, Mie, Jue, Vie);
        } else {
            iss >> P_r; iss >> P_p;iss >> P_u;
            // cout << P_r <<' '<<P_p <<' '<<P_u << endl;
        }
        c++;
    }
    M_maquinas_o = M_maquinas; N_pacientes_o= N_pacientes; D_doctores_o=D_doctores; 
    P_r_o=P_r; P_p_o= P_p; P_u_o=P_u;
    archivo.close();
}


void DList(vector<int> doctors_schedule[5], int Lun, int Ma, int Mie, int Jue, int Vie) {
    doctors_schedule[0].push_back(Lun);
    doctors_schedule[1].push_back(Ma);
    doctors_schedule[2].push_back(Mie);
    doctors_schedule[3].push_back(Jue);
    doctors_schedule[4].push_back(Vie);
}




void movement(vector<int[][20*16]> neigborhood, int S_c[][20*16],queue<int> &Tabu_list, int U, int P, int R) {
        for (size_t i = 0; i < U+P+R; i++)
        {
            if (i in list(Tabu_list.queue)) {continue;}

            // S_c_aux = deepcopy(S_c) CHECAR ESTO
            int S_c_aux[U+P+R][20*16];
            int aux[20*16];

            S_c_aux = S_c; // hace un deep copy
            aux = S_c_aux[i];
            if(i == U-1) {
                S_c_aux[i] = S_c_aux[0];
                S_c_aux[0] = aux;
                neigborhood.push_back(S_c_aux)
            } else if(i == U+P-1) {
                S_c_aux[i] = S_c_aux[U];
                S_c_aux[U] = aux;
                neigborhood.push_back(S_c_aux)
            } else if(i == U+P+R-1) {
                S_c_aux[i] = S_c_aux[P+U];
                S_c_aux[P+U] = aux;
                neigborhood.push_back(S_c_aux)
            } else {
                S_c_aux[i] = S_c_aux[i+1];
                S_c_aux[i+1] = aux;
                neigborhood.push_back(S_c_aux);
            }
        }
    }

int f_ev(int sol_d[][20*16], int sol_m[][20*16], int N_pac) {
        int same_doctor = 0;
        int same_machine = 0;
        for (size_t index_pac = 0; index_pac < N_pac; index_pac++) { // sol_d
            vector<int> doctores_pac;
            vector<int> machines_pac;
            for (size_t index_doc = 0; index_pac < N_pac; index_pac++) { // sol_d[index_pac]
                if (index_doc != 0) {
                    if (doctor in doctores_pac) { same_doctor = same_doctor +1};
                    if (sol_m[index_pac][index_doc] in machines_pac) { same_machine = same_machine +1};
                    doctores_pac.push_back(doctor);
                    machines_pac.push_back(sol_m[index_pac][index_doc]);
                    continue //#in C++ skip to another daylight
                }
            }
        }
        return same_doctor+same_machine
}

void best_N(int S_c[][20*16], int S_m[][20*16], queue<int> &Tabu_list, int U, int P, int R) {
    int local_best_d[][20*16];
    int local_best_m[][20*16];
    vector<int[][20*16]> neigborhood_d;
    vector<int[][20*16]> neigborhood_m;
    movement(neigborhood_d, S_c, Tabu_list, U,P,R);
    movement(neigborhood_m, S_c, Tabu_list, U,P,R);

    local_best_d = neigborhood_d[0];
    local_best_m = neigborhood_m[0];
    int forbidden_movement = 0;
    for (size_t index = 0; index < neigborhood_d.size(); index++){
        if (f_ev(neigborhood_d[index], neigborhood_m[index], U+P+R) > f_ev(local_best_d, neigborhood_m[index], U+P+R)){
            local_best_d = neigborhood_d[index];
            local_best_m = neigborhood_m[index];
            forbidden_movement = index;
        }
    }
    return local_best_d, local_best_m, forbidden_movement
}

void update_tabu_list (queue<int> &Tabu_list, int forbidden_movement) {
        if(not Tabu_list.full(){
            Tabu_list.push(forbidden_movement);
        } else {
            cout << " Tabu list front elemnt " << Tabu_list.front() << endl;
            Q_paciente.pop();
            Tabu_list.push(forbidden_movement);
        }
}

void Tabu(int S_c[][20*16], int S_m[][20*16], int max_size, int U, int P, int R, int stop_criterion) {
    queue<int> Tabu_list; // DEFINIR MAX SIZE
    int S_best_d[U+P+R][20*16]; int S_best_m[U+P+R][20*16];
    int sd[U+P+R][20*16]; int sm[U+P+R][20*16];
    int S_v_d[U+P+R][20*16]; int S_v_m[U+P+R][20*16];
    int forbidden_movement;

    S_best_d = S_c_d;
    S_best_m = S_c_m;

    sd, sm,_= best_N(S_c_d, S_c_m, Tabu_list, U,P,R);
    int c=0;
    printf("stop_criterion: ", stop_criterion, "max size: ", max_size);
    while (c!=stop_criterion) {
        S_v_d, S_v_m, forbidden_movement = best_N(S_c_d, S_c_m, Tabu_list, U,P,R);
        S_c_d = S_v_d;
        S_c_m = S_v_m;
        update_tabu_list(Tabu_list, forbidden_movement);
        if (f_ev(S_c_d, S_c_m, U,P,R) > f_ev(S_best_d, S_best_m, U,P,R)) {
            S_best_d = S_c_d;
            S_best_m = S_c_m;
        }
        c +=1;
    }
    printf("c: ", c);
    return S_best_d, S_best_m;
}
