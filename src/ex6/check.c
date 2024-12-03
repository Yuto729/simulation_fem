#include <stdio.h>
#include <stdlib.h>

// 力の成分を計算する関数
void calculateForceSum(const char *filename) {
    FILE *file;
    char buffer[256];
    double fx_sum = 0.0, fy_sum = 0.0, fz_sum = 0.0;
    double fx, fy, fz;

    // ファイルを開く
    file = fopen(filename, "r");
    if (file == NULL) {
        printf("Error: Cannot open file %s\n", filename);
        return;
    }

    // ヘッダー行をスキップ
    fgets(buffer, sizeof(buffer), file);

    // 各行を読み取り、力の総和を計算
    while (fgets(buffer, sizeof(buffer), file)) {
        // ファイルのフォーマットに従い、力成分を読み取る
        sscanf(buffer, "%*d,%*lf,%*lf,%*lf,%*lf,%*lf,%*lf,%lf,%lf,%lf", &fx, &fy, &fz);
        fx_sum += fx;
        fy_sum += fy;
        fz_sum += fz;
    }

    // ファイルを閉じる
    fclose(file);

    // 結果を出力
    printf("Total Force Sum:\n");
    printf("Fx: %.10e\n", fx_sum);
    printf("Fy: %.10e\n", fy_sum);
    printf("Fz: %.10e\n", fz_sum);

    // 釣り合いの状態を確認
    if (fabs(fx_sum) < 1e-6 && fabs(fy_sum) < 1e-6 && fabs(fz_sum) < 1e-6) {
        printf("The forces are in equilibrium.\n");
    } else {
        printf("The forces are not in equilibrium.\n");
    }
}

int main(int argc, char *argv[]) {

    calculateForceSum("result.csv");
    return 0;
}
