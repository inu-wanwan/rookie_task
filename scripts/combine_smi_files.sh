#!/bin/bash

# 出力ファイル名
output_file="gen_all.smi"

# 既存の出力ファイルを削除（存在する場合）
if [ -f "$output_file" ]; then
    rm "$output_file"
fi

# 全ての入力ファイルを結合
for smi_file in "$@"
do
    cat "$smi_file" >> "$output_file"
    echo "" >> "$output_file"  # ファイル間に空行を追加
done

echo "全てのSMILESファイルが $output_file にまとめられました。"
