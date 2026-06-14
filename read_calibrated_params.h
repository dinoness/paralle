/**
 * read_calibrated_params.h — C 语言读取 calibrated_params.csv 的规则与示例
 *
 * CSV 格式说明：
 *   - 以 '#' 开头的行为注释，跳过
 *   - 数据行格式: "标签,值1[,值2,值3]"
 *   - 标签前缀规则:
 *       r1, r2, alpha       → 单值（标量）
 *       l0_N                → 单值，N ∈ [1,5]
 *       B_N, Pm_N           → 三值 (x,y,z)，N ∈ [1,5]
 *       limb_dir_N          → 二值 (theta_base, theta_move)，N ∈ [1,5]
 *   - 所有值使用 SI 单位: 长度 [m], 角度 [rad]
 *   - 精度: 12 位小数，C 中 double 可完整保存
 *
 * 使用方法:
 *   1. 将 calibrated_params.csv 放在可执行文件同目录
 *   2. 调用 read_calibrated_params() 加载全部参数
 */

#ifndef READ_CALIBRATED_PARAMS_H
#define READ_CALIBRATED_PARAMS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* 机构参数结构体 */
typedef struct {
    double r1;                     /* 动平台下半径 [m] */
    double r2;                     /* 动平台上半径 [m] */
    double alpha;                  /* U 副倾角 [rad] */
    double l0[5];                  /* 初始杆长 [m] */
    double B[5][3];                /* 基座铰点 [m] */
    double Pm[5][3];               /* 动平台铰点 (平台系) [m] */
    double limb_dir[5][2];         /* 支链方向角 [θ_base, θ_move] [rad] */
} CalibratedParams;

/**
 * 从 calibrated_params.csv 加载标定参数
 * 返回 0 成功，-1 失败
 */
static int read_calibrated_params(const char *filename, CalibratedParams *p) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "read_calibrated_params: cannot open %s\n", filename);
        return -1;
    }

    char line[256];
    int l0_idx = 0, B_idx = 0, Pm_idx = 0, limb_idx = 0;

    while (fgets(line, sizeof(line), fp)) {
        /* 跳过注释和空行 */
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;

        /* 去除末尾换行符 */
        size_t len = strlen(line);
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r'))
            line[--len] = '\0';

        char label[32];
        double v1 = 0, v2 = 0, v3 = 0;
        int n = sscanf(line, "%31[^,],%lf,%lf,%lf", label, &v1, &v2, &v3);

        if (n < 2) continue;  /* 格式错误 */

        /* ── 标量 ── */
        if      (strcmp(label, "r1")    == 0) p->r1    = v1;
        else if (strcmp(label, "r2")    == 0) p->r2    = v1;
        else if (strcmp(label, "alpha") == 0) p->alpha = v1;

        /* ── l0_N ── */
        else if (strncmp(label, "l0_", 3) == 0) {
            int idx = atoi(label + 3);
            if (idx >= 1 && idx <= 5) { p->l0[idx-1] = v1; l0_idx++; }
        }

        /* ── B_N (x,y,z) ── */
        else if (strncmp(label, "B_", 2) == 0 && !strchr(label, 'm')) {
            int idx = atoi(label + 2);
            if (idx >= 1 && idx <= 5 && n >= 4) {
                p->B[idx-1][0] = v1; p->B[idx-1][1] = v2; p->B[idx-1][2] = v3;
                B_idx++;
            }
        }

        /* ── Pm_N (x,y,z) ── */
        else if (strncmp(label, "Pm_", 3) == 0) {
            int idx = atoi(label + 3);
            if (idx >= 1 && idx <= 5 && n >= 4) {
                p->Pm[idx-1][0] = v1; p->Pm[idx-1][1] = v2; p->Pm[idx-1][2] = v3;
                Pm_idx++;
            }
        }

        /* ── limb_dir_N (θ_base, θ_move) ── */
        else if (strncmp(label, "limb_dir_", 9) == 0) {
            int idx = atoi(label + 9);
            if (idx >= 1 && idx <= 5 && n >= 3) {
                p->limb_dir[idx-1][0] = v1;
                p->limb_dir[idx-1][1] = v2;
                limb_idx++;
            }
        }
    }

    fclose(fp);

    /* 完整性校验 */
    if (l0_idx != 5 || B_idx != 5 || Pm_idx != 5 || limb_idx != 5) {
        fprintf(stderr, "read_calibrated_params: incomplete data "
                "(l0=%d/5 B=%d/5 Pm=%d/5 limb=%d/5)\n",
                l0_idx, B_idx, Pm_idx, limb_idx);
        return -1;
    }

    return 0;
}

/*
 * ── 使用示例 ──
 *
 * int main(void) {
 *     CalibratedParams p;
 *     if (read_calibrated_params("calibrated_params.csv", &p) != 0)
 *         return 1;
 *
 *     printf("r1 = %.6f mm\n",  p.r1 * 1000.0);
 *     printf("r2 = %.6f mm\n",  p.r2 * 1000.0);
 *     printf("alpha = %.4f deg\n", p.alpha * 180.0 / M_PI);
 *     for (int i = 0; i < 5; i++) {
 *         printf("limb %d: l0=%.4f mm  B=(%.3f,%.3f,%.3f) mm  "
 *                "Pm=(%.3f,%.3f,%.3f) mm  "
 *                "dir=(%.2f, %.2f) deg\n",
 *                i+1, p.l0[i]*1000,
 *                p.B[i][0]*1000, p.B[i][1]*1000, p.B[i][2]*1000,
 *                p.Pm[i][0]*1000, p.Pm[i][1]*1000, p.Pm[i][2]*1000,
 *                p.limb_dir[i][0]*180/M_PI, p.limb_dir[i][1]*180/M_PI);
 *     }
 *     return 0;
 * }
 */

#endif /* READ_CALIBRATED_PARAMS_H */

