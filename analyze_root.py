#!/usr/bin/env python3

import os
import sys

def analyze_root_file(filepath):
    """루트 파일의 기본 정보를 분석합니다."""
    if not os.path.exists(filepath):
        print(f"파일이 존재하지 않습니다: {filepath}")
        return
    
    file_size = os.path.getsize(filepath)
    print(f"파일: {filepath}")
    print(f"크기: {file_size} bytes")
    
    if file_size == 0:
        print("경고: 파일이 비어있습니다!")
    elif file_size < 1000:
        print("경고: 파일이 매우 작습니다!")
    else:
        print("파일 크기가 정상입니다.")

if __name__ == "__main__":
    # 기존 파일과 새 파일 비교
    old_file = "/u/user/rmsvlf000/Trigger/test_0.root"
    new_file = "/u/user/rmsvlf000/Trigger/build/test_fixed_0.root"
    
    print("=== 루트 파일 분석 ===")
    print()
    
    print("1. 기존 파일 (마이그레이션 전):")
    analyze_root_file(old_file)
    print()
    
    print("2. 수정된 파일 (마이그레이션 후):")
    analyze_root_file(new_file)
    print()
    
    # 크기 비교
    if os.path.exists(old_file) and os.path.exists(new_file):
        old_size = os.path.getsize(old_file)
        new_size = os.path.getsize(new_file)
        
        print("=== 크기 비교 ===")
        print(f"기존 파일: {old_size} bytes")
        print(f"새 파일: {new_size} bytes")
        
        if new_size > old_size:
            print("✅ 새 파일이 더 큽니다 - 데이터가 더 많이 저장되었습니다!")
        elif new_size == old_size:
            print("⚠️  파일 크기가 동일합니다")
        else:
            print("❌ 새 파일이 더 작습니다 - 데이터 손실 가능성")

