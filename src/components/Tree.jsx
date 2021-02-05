import SortableTree from 'react-sortable-tree'


const Tree = props => (
  <ul>
    {props.rules.map(rule => (
      <li key={rule.uuid}>
        {rule.name}
        {rule.children ? <Tree rules={rule.children}/> : null}
      </li>
    ))}
  </ul>
)

export const RuleTree = props => {
  return (
    <SortableTree
      treeData={props.rules}
      isVirtualized={false}
      onChange={props.onChange}
    />
  )
}
