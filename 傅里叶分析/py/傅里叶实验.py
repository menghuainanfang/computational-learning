"""第24周两组傅里叶实验：控制台、CSV、JSON 与时域/频域对比图。"""
from pathlib import Path
from time import perf_counter
import argparse
import csv
import json
import platform
import sys
from datetime import datetime
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def _序列(x):
    a = np.asarray(x, dtype=complex)
    if a.ndim != 1 or a.size == 0:
        raise ValueError("输入必须是非空一维序列")
    return a


def dft(x):
    """直接按定义逐频率求和；正变换不除以 N。"""
    a = _序列(x)
    n = np.arange(a.size)
    return np.array([np.sum(a * np.exp(-2j * np.pi * k * n / a.size))
                     for k in range(a.size)])


def fft_radix2(x):
    """递归拆分偶下标与奇下标，合并为同一个 DFT。"""
    a = _序列(x)
    if a.size & (a.size - 1):
        raise ValueError("基-2 FFT 的长度必须为正的 2 的幂")

    def 递归(v):
        n = v.size
        if n == 1:
            return v.copy()
        偶部 = 递归(v[::2])
        奇部 = 递归(v[1::2])
        旋转奇部 = np.exp(-2j * np.pi * np.arange(n // 2) / n) * 奇部
        return np.concatenate((偶部 + 旋转奇部, 偶部 - 旋转奇部))

    return 递归(a)


def _计时(方法, x, 次数):
    # 预热不计入；计时内仅执行变换，包含算法自身的分配和输入检查。
    方法(x)
    开始 = perf_counter()
    for _ in range(次数):
        y = 方法(x)
    return y, (perf_counter() - 开始) / 次数


def _写表(路径, 表):
    with 路径.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(表[0]))
        w.writeheader()
        w.writerows(表)


def main():
    parser = argparse.ArgumentParser(description="运行第24周两组傅里叶实验")
    parser.add_argument("--no-plots", action="store_true", help="只输出数值数据，不生成图像")
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[1] / "results",
                        help="结果输出目录，默认个人目录下的results")
    args = parser.parse_args()
    输出 = args.output
    输出.mkdir(parents=True, exist_ok=True)
    容差 = 1e-10
    次数 = 10
    随机源 = np.random.default_rng(20261109)
    变换表 = []
    for n in (8, 64, 256):
        x = 随机源.standard_normal(n) + 1j * 随机源.standard_normal(n)
        参考 = np.fft.fft(x)
        for 名称, 方法 in (("直接DFT", dft), ("基2FFT", fft_radix2)):
            y, 时间 = _计时(方法, x, 次数)
            误差 = float(np.max(np.abs(y - 参考)))
            assert 误差 < 容差, (名称, n, 误差)
            变换表.append({"N": n, "方法": 名称, "最大绝对误差": 误差,
                          "平均耗时_s": 时间, "重复次数": 次数})

    # 实序列的负频率下标由 (-k) % N 得到，直流与 Nyquist 自共轭。
    实序列 = np.array([1., 2., -1., 0., 3., -2., 4., 1.])
    频谱 = fft_radix2(实序列)
    对称误差 = float(np.max(np.abs(频谱[(-np.arange(8)) % 8] - 频谱.conj())))
    时域能量 = float(np.sum(实序列 ** 2))
    频域能量 = float(np.sum(np.abs(频谱) ** 2) / 8)
    重构误差 = float(np.max(np.abs(np.fft.ifft(频谱) - 实序列)))
    性质 = {"共轭对称误差": 对称误差, "时域能量": 时域能量,
            "频域能量": 频域能量, "能量差": abs(时域能量 - 频域能量),
            "重构误差": 重构误差}
    assert max(对称误差, 性质["能量差"], 重构误差) < 容差
    # 另外用同一实序列核对讲义中的线性及循环右移相位约定。
    z = np.arange(8.)
    性质["线性误差"] = float(np.max(np.abs(
        fft_radix2(2 * 实序列 - 3 * z) - (2 * 频谱 - 3 * fft_radix2(z)))))
    性质["循环右移误差"] = float(np.max(np.abs(
        fft_radix2(np.roll(实序列, 2)) - np.exp(-2j * np.pi * np.arange(8) * 2 / 8) * 频谱)))
    assert max(性质["线性误差"], 性质["循环右移误差"]) < 容差
    for 非法 in ([], [1, 2, 3], np.ones((2, 2))):
        try:
            fft_radix2(非法)
        except ValueError:
            pass
        else:
            raise AssertionError("基2 FFT 未拒绝非法长度或维度")

    fs, n = 128., 256
    t = np.arange(n) / fs
    目标 = np.sin(2 * np.pi * 5 * t)
    原信号 = 目标 + .3 * np.sin(2 * np.pi * 30 * t)
    频率 = np.fft.rfftfreq(n, d=1 / fs)
    X = np.fft.rfft(原信号)
    幅值 = np.abs(X) / n
    幅值[1:-1] *= 2  # 本例 N 为偶数，末项为 Nyquist，不翻倍。
    主峰 = np.argsort(幅值[1:])[-2:] + 1
    主峰 = 主峰[np.argsort(频率[主峰])]
    assert np.allclose(频率[主峰], [5, 30], rtol=0, atol=容差)
    assert np.allclose(幅值[主峰], [1, .3], rtol=0, atol=容差)
    Y = X.copy()
    Y[频率 > 10] = 0
    滤波 = np.fft.irfft(Y, n=n)
    滤后幅值 = np.abs(Y) / n
    滤后幅值[1:-1] *= 2
    rmse前 = float(np.sqrt(np.mean((原信号 - 目标) ** 2)))
    rmse后 = float(np.sqrt(np.mean((滤波 - 目标) ** 2)))
    assert abs(rmse前 - .3 / np.sqrt(2)) < 容差
    assert rmse后 < 容差
    时域表 = [{"n": i, "时间_s": float(t[i]), "目标": float(目标[i]),
               "滤波前": float(原信号[i]), "滤波后": float(滤波[i])} for i in range(n)]
    频域表 = [{"频率_Hz": float(f), "滤波前单边幅值": float(a),
               "滤波后单边幅值": float(b)} for f, a, b in zip(频率, 幅值, 滤后幅值)]
    a, b = [1, 2, 3], [1, 1]
    手算 = np.array([1., 3., 5., 3.])
    卷积 = np.fft.ifft(np.fft.fft(a, 4) * np.fft.fft(b, 4))
    卷积误差 = float(np.max(np.abs(卷积 - 手算)))
    库误差 = float(np.max(np.abs(卷积 - np.convolve(a, b))))
    assert max(卷积误差, 库误差) < 容差
    if not args.no_plots:
        # 三幅图分别展示波形、恢复误差和单边幅值，均使用本次实验的数组。
        fig, axes = plt.subplots(3, 1, figsize=(10, 9), constrained_layout=True)
        axes[0].plot(t, 原信号, label="Input", alpha=.6)
        axes[0].plot(t, 目标, label="Target 5 Hz", color="black", linestyle="--")
        axes[0].plot(t, 滤波, label="Low-pass output", linewidth=1)
        axes[0].set(xlabel="Time (s)", ylabel="Amplitude", xlim=(0, .5))
        axes[1].plot(t, 原信号-目标, label="Before - target")
        axes[1].plot(t, 滤波-目标, label="After - target")
        axes[1].set(xlabel="Time (s)", ylabel="Error", xlim=(0, .5))
        axes[2].plot(频率, 幅值, label="Before", marker=".")
        axes[2].plot(频率, 滤后幅值, label="After", linestyle="--")
        axes[2].axvline(10, color="gray", linestyle=":", label="Cutoff 10 Hz")
        axes[2].set(xlabel="Frequency (Hz)", ylabel="One-sided amplitude", xlim=(0,64))
        for ax in axes:
            ax.grid(alpha=.2)
            ax.legend()
        fig.savefig(输出 / "实验B_对比.svg")
        fig.savefig(输出 / "实验B_对比.png", dpi=160)
        plt.close(fig)
    汇总 = {"运行日期": datetime.now().astimezone().isoformat(),
            "Python可执行文件": sys.executable, "Matplotlib": matplotlib.__version__,
            "Python": platform.python_version(), "NumPy": np.__version__,
            "随机种子": 20261109, "绝对容差": 容差,
            "实验A": 变换表, "性质验证": 性质,
            "实验B": {"fs_Hz": fs, "N": n, "频率间隔_Hz": fs / n,
                      "主要频率_Hz": 频率[主峰].tolist(), "单边幅值": 幅值[主峰].tolist(),
                      "滤波前RMSE": rmse前, "滤波后RMSE": rmse后,
                      "FFT线性卷积": 卷积.real.tolist(), "相对手算最大误差": 卷积误差,
                      "相对numpy卷积最大误差": 库误差}}
    _写表(输出 / "实验A_变换误差与计时.csv", 变换表)
    _写表(输出 / "实验B_时域对照.csv", 时域表)
    _写表(输出 / "实验B_频域对照.csv", 频域表)
    (输出 / "实验结果.json").write_text(json.dumps(汇总, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(汇总, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
