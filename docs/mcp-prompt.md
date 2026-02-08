# MCP Tool Implementation Prompt

## Репозиторий
`d:\Prototype\SuffixTree`

## Чеклист
`docs/mcp-checklist.md` — секция 4 (Per-Tool Checklist). Это единственный источник правды о том, что реализовано и что нет.

## Эталонные файлы
- Tool + records: `src/Seqeron/Mcp/Seqeron.Mcp.Parsers/Tools/ParsersTools.cs`
- Тесты: `src/Seqeron/Mcp/Seqeron.Mcp.Parsers.Tests/`
- Документация: `docs/mcp/tools/parsers/*.md`, `*.mcp.json`

## Цикл на каждый tool

1. Найди следующий невыполненный tool в `docs/mcp-checklist.md`
2. Изучи существующую реализацию в `Seqeron.Genomics` (DocRef из чеклиста)
3. Реализуй tool + result record по образцу эталонных файлов
4. Создай тесты (Schema + Binding)
5. Создай документацию (.md + .mcp.json)
6. `dotnet build && dotnet test` — все тесты должны пройти
7. Отметь `[x]` пункты a-j в `docs/mcp-checklist.md`
8. Коммит: `feat(MCP/Parsers): add <tool_name> tool`
9. Повтори с шага 1

## Критически важно

- Mission-critical библиотека. Ошибки и неточности недопустимы.
- Каждый tool должен корректно вызывать метод из `Seqeron.Genomics`, а не реализовывать логику заново.
- Тесты должны проверять реальное поведение, а не заглушки.

## Задача

Начни со следующего невыполненного tool из чеклиста.
